"""
Tools for managing metrics with metadata attached.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

from __future__ import division
import pandas as pd
import numpy as np
import os
import martian

# Constants for computing CTCF footprinting score
UPSTR_CTCF_PEAK_POS = '-23'
DOWNSTR_CTCF_PEAK_POS = '34'
CTCF_VALLEY_POS = '-3'

class MetricAnnotations:
    def __init__(self):
        """Load in metric information from the associated csv file.
        """
        file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'metrics.csv')
        data = pd.read_csv(file_path, index_col=0)
        # Makes empty values None instead of NaN
        data = data.where(pd.notnull(data), None)
        self.metric_data = {}
        for key in data.index.values:
            self.metric_data[key] = data.loc[key]

    def gen_metric(self, key, value, species=None, debug=True):
        """Returns a single metric object for the given key and value.  Alerts are raised when metric falls outside
        normal range, which depends on debug status."""
        metric_info = self.metric_data[key]
        name = metric_info.full_name
        alert_name = metric_info.alert_name
        if species is not None and species:
            key += '_{species}'.format(**locals())
            name += ' ({species})'.format(**locals())
            # alert_name may not be specified for some metrics, in which case an alert can't be raised, therefore
            # no need to worry about its species specific name format.
            if alert_name is not None:
                alert_name += ' ({species})'.format(**locals())

        # alarm ranges are dependent on debug, which indicates internal or customer-facing.
        description = (metric_info.description if debug else metric_info.description_cs)
        acceptable = (metric_info.acceptable if debug else metric_info.acceptable_cs)
        targeted = (metric_info.targeted if debug else metric_info.targeted_cs)

        kwargs = {
            'description': description,
            'acceptable': acceptable,
            'targeted': targeted,
            'evaluate_type': metric_info.evaluate_type,
            'format_type': metric_info.format_type,
            'category': metric_info.category,
            'alert_name': alert_name,
        }
        return Metric(key, name, value, **kwargs)

    def compile_summary_metrics(self, value_dict, keys=None, species_list=None):
        """Processes a metrics dictionary and select summary metrics based on 2nd column in metrics.csv
        for keys provided or all registered keys in metrics.csv"""
        keylist = self.metric_data.keys() if keys is None else keys
        output = {}
        for key in keylist:
            if key in self.metric_data:
                metric_info = self.metric_data[key]
                # if it is a summary metric
                if metric_info.summary:
                    # do species specific processing
                    if metric_info.is_species_specific:
                        if species_list is None or not species_list:
                            raise ValueError('Must provide a species list for species-specific metrics')
                        for species in species_list:
                            key_suffix = "" if len(species_list) == 1 else "_{}".format(species)
                            subkey = '{key}{key_suffix}'.format(**locals())
                            if subkey in value_dict:
                                output.update({subkey: value_dict[subkey]})
                            else:
                                martian.log_info('{} not found in metrics'.format(key))
                    else:
                        if key in value_dict:
                            output.update({key: value_dict[key]})
                        else:
                            martian.log_info('{} not found in metrics'.format(key))
            else:
                martian.log_info('{} not found in registered metrics'.format(key))
        return output

    def gen_metric_helptext(self, keys):
        """Processes a metrics dictionary and generates helptext for keys if present in metrics.csv"""
        output = []
        for key in keys:
            if key in self.metric_data:
                metric_info = self.metric_data[key]
                if metric_info.help_description is not None:
                    output += [[metric_info.full_name, [metric_info.help_description]]]
            else:
                martian.log_info('{} not found in registered metrics'.format(key))
        return output

    def gen_metric_list(self, value_dict, keys, species_list=None, species_is_library=False, debug=True):
        """Returns a list of metric objects for the provided keys, using the value dictionary to get values for
        each metric.  When metrics are species-specific, a list of species is required.  Alerts are raised when
        metrics fall outside normal ranges, which depend on debug status."""
        output = []
        for key in keys:
            metric_info = self.metric_data[key]
            if metric_info.is_species_specific:
                if species_list is None or not species_list:
                    raise ValueError('Must provide a species list for species-specific metrics')
                for species in species_list:
                    key_suffix = "" if len(species_list) == 1 else "_{}".format(species)
                    if species_is_library:
                        key_suffix = "_{}".format(species)
                    subkey = '{key}{key_suffix}'.format(**locals())
                    if subkey in value_dict:
                        output.append(self.gen_metric(key, value_dict[subkey], species, debug))
                    else:
                        martian.log_info('{} not found in metrics'.format(subkey))
            else:
                if key in value_dict:
                    output.append(self.gen_metric(key, value_dict[key], debug=debug))
                else:
                    martian.log_info('{} not found in metrics'.format(key))
        return output


class Metric:
    """Contains metadata about a single metric along with methods for evaluating its value with respect to that
    metadata.
    """
    def __init__(self, key, name, value,
                 description=None, acceptable=None, targeted=None, evaluate_type=None,
                 format_type=None, category=None, alert_name=None):
        self.key = key
        self.name = name
        self.value = value
        self.description = description
        try:
            self.acceptable = float(acceptable)
        except (ValueError, TypeError):
            self.acceptable = acceptable
        try:
            self.targeted = float(targeted)
        except (ValueError, TypeError):
            self.targeted = targeted
        self.alert_name = alert_name

        if evaluate_type not in [None, 'lt', 'gt', 'range', 'exists', 'is_true', 'is_false']:
            raise ValueError('Unknown evaluation type: {}'.format(evaluate_type))
        function_map = {None: lambda x: True,
                        'lt': self._less_than,
                        'gt': self._greater_than,
                        'range': self._in_range,
                        'exists': self._exists,
                        'is_true': self._is_true,
                        'is_false': self._is_false}
        self.evaluate_type = evaluate_type
        self.evaluation_function = function_map[evaluate_type]

        if format_type is None:
            format_type = 'flat'
        if format_type not in ['flat', 'float', 'percentage', 'int']:
            raise ValueError('Unknown format type: {}'.format(format_type))
        self.format_type = format_type

        if category is None:
            category = 'General'
        self.category = category

    def gen_metric_dict(self, default_threshold=None):
        if default_threshold is not None:
            assert default_threshold in ['pass', 'warn', 'error']
        threshold = self.threshold_type
        translate_dict = {
            'VALID': 'pass',
            'WARN': 'warn',
            'ERROR': 'error',
        }
        return {
            "threshold": translate_dict[threshold] if default_threshold is None else default_threshold,
            "metric": self.value_string,
            "name": self.name,
        }

    @property
    def alarm_dict(self):
        threshold = self.threshold_type
        if threshold == "VALID":
            return {}

        return {
            "raw_value": self.value,
            "formatted_value": self.value_string,
            "raised": True,
            "parent": self.key,
            "title": self.alert_name,
            "message": self.description,
            "level": threshold,
            "test": "",
            "id": self.key
        }

    # Evaluation methods
    @staticmethod
    def _less_than(value, target):
        return value < target

    @staticmethod
    def _greater_than(value, target):
        return value > target

    @staticmethod
    def _in_range(value, target):
        return (value > target[0]) and (value < target[1])

    @staticmethod
    def _exists(value, target):
        return value is not None

    @staticmethod
    def _is_true(value, target):
        return value is True

    @staticmethod
    def _is_false(value, target):
        return value is False

    @property
    def threshold_type(self):
        if self.targeted is None and self.acceptable is None:
            return "VALID"
        elif self.targeted is None:
            # Only acceptable - error if we don't meet them.
            if self.evaluation_function(self.value, self.acceptable):
                return "VALID"
            return "ERROR"
        elif self.acceptable is None:
            # Only targets - warn if we don't meet them.
            if self.evaluation_function(self.value, self.targeted):
                return "VALID"
            return "WARN"
        else:
            # Both set - error/warn depending on which we meet.
            if self.evaluation_function(self.value, self.targeted):
                return "VALID"
            elif self.evaluation_function(self.value, self.acceptable):
                return "WARN"
            else:
                return "ERROR"

    @property
    def color(self):
        if self.targeted is None and self.acceptable is None:
            # Meaningless, return grey hex code
            return 'BEBEBE'
        threshold = self.threshold_type
        if threshold == "VALID":
            return "B4FFB4"
        elif threshold == "WARN":
            return "FFFFB4"
        else:
            return "FFB4B4"

    @property
    def acceptable_string(self):
        return self._format_target_value(self.acceptable)

    @property
    def targeted_string(self):
        return self._format_target_value(self.targeted)

    @property
    def value_string(self):
        return self._format_value(self.value)

    def _format_target_value(self, value):
        if value is None:
            return ''

        if self.evaluate_type == "exists":
            return 'Exists'
        if self.evaluate_type == "is_true":
            return 'True'
        if self.evaluate_type == "is_false":
            return 'False'

        if self.evaluate_type == "lt":
            return '< {}'.format(self._format_value(value))
        if self.evaluate_type == "gt":
            return '> {}'.format(self._format_value(value))
        if self.evaluate_type == "range":
            return '{} - {}'.format(self._format_value(value[0]), self._format_value(value[1]))

    def _format_value(self, value):
        if value is None:
            return 'None'

        if self.format_type == 'flat':
            return '{}'.format(value)
        elif self.format_type == 'float':
            return '{:,.2f}'.format(value)
        elif self.format_type == 'percentage':
            return '{:.1%}'.format(value)
        elif self.format_type == 'int':
            return '{:,.0f}'.format(value)


def calculate_tss_score_and_profile(relative_positions):
    tss_df = pd.read_csv(relative_positions)
    xvals = np.arange(-1000, 1000 + 1)
    yvals = tss_df[map(str, xvals)].sum(axis=0)
    if any(yvals > 0):
        min_nonzero = yvals[yvals > 0].min()
        yvals /= min_nonzero
        score = yvals.max()
        return score, np.array(yvals), xvals
    else:
        return 0.0, np.array(yvals), xvals


def calculate_ctcf_score_and_profile(relative_positions):
    ctcf_df = pd.read_csv(relative_positions)
    xvals = np.arange(-250, 250 + 1)
    yvals = ctcf_df[map(str, xvals)].sum(axis=0)

    if max(yvals) > 0:
        yvals /= max(yvals)
        score = yvals[UPSTR_CTCF_PEAK_POS] + yvals[DOWNSTR_CTCF_PEAK_POS] - yvals[CTCF_VALLEY_POS]
        return score, np.array(yvals), xvals
    else:
        return 0.0, np.array(yvals), xvals

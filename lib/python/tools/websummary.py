"""
Tools for quickly creating a simple websummary, with embedded figures.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

import os
import base64
from StringIO import StringIO

CSS = """
/* Style the metric tables */
table {
    border-collapse: collapse;
}
th, td {
    padding: 15px;
}
tr:nth-child(even) {background-color: #f2f2f2}
table, th, td {
    border: 1px solid black;
}

/* helvetica-style all the text */
h1, h3, p, blockquote, tr, pre, table, th, td, ul, li{
    font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
    font-style: normal;
    font-variant: normal;
}

h1 {
    text-align: center;
    font-size: 24px;
    font-weight: 500;
    line-height: 26.4px;
}

p, blockquote{
    font-size: 14px;
    font-weight: 400;
    line-height: 20px;
}

.matplotlib, .metrics{
    text-align: center;
}

/* Style the buttons that are used to open and close the accordion panel */
.accordion {
    background-color: #eee;
    color: #444;
    cursor: pointer;
    padding: 18px;
    width: 100%;
    font-size: large;
    text-align: left;
    border-style: solid;
    border-width: 1px;
    outline: none;
    transition: 0.4s;
}

/* Add a background color to the button if it is clicked on (add the .active class with JS), 
and when you move the mouse over it (hover) */
.active, .accordion:hover {
    background-color: #ccc;
}

/* Style the accordion panel. Note: hidden by default */
.panel {
    padding: 0 18px;
    border-style: solid;
    border-width: 1px;
    background-color: white;
    display: none;
}

/* Style the tab */
.tab {
    overflow: hidden;
    border: 1px solid #ccc;
    background-color: #cce;
}

/* Style the buttons that are used to open the tab content */
.tab button {
    font-size: large;
    background-color: #cce;
    float: left;
    border: none;
    outline: none;
    cursor: pointer;
    padding: 14px 16px;
    transition: 0.3s;
}

/* Change background color of buttons on hover */
.tab button:hover {
    background-color: #ddf;
}

/* Create an active/current tablink class */
.tab button.active {
    background-color: #bbd;
}

/* Style the tab content */
.tabcontent {
    display: none;
    padding: 0px;
    border: none;
}
"""

SCRIPT = """
var acc = document.getElementsByClassName("accordion");
var i;

for (i = 0; i < acc.length; i++) {
    acc[i].onclick = function(){
        this.classList.toggle("active");
        var panel = this.nextElementSibling;
        if (panel.style.display === "block") {
            panel.style.display = "none";
        } else {
            panel.style.display = "block";
        }
    }
}

function openTab(evt, tabName) {
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
    }

    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }

    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.className += " active";
}

// Make sure the default element is open on load
document.getElementById("defaultOpen").click();
"""


def metric_list_to_html(metric_list):
    data = '\n'.join('<tr>'
                     '<td align="right">{name}</td>'
                     '<td align="right">{value}</td>'
                     '</tr>'.format(name=metric.name, value=metric.value_string)
                     for metric in metric_list)
    return """
    <table border="1" class="metric_table">
    {}
    </table>
    """.format(data)


def targeted_metric_list_to_html(metric_list):
    data = '\n'.join('<tr style="background-color:#{color}">'
                     '<td align="right">{category}</td>'
                     '<td align="right">{name}</td>'
                     '<td align="right">{value}</td>'
                     '<td align="right">{acceptable}</td>'
                     '<td align="right">{targeted}</td>'
                     '</tr>'.format(color=metric.color, category=metric.category, name=metric.name,
                                    value=metric.value_string, acceptable=metric.acceptable_string,
                                    targeted=metric.targeted_string)
                     for metric in metric_list)
    return """
    <table border="1" class="metric_table">
    <th align="right">Category</th>
    <th align="right">Metric</th>
    <th align="right">Observed Value</th>
    <th align="right">Acceptable</th>
    <th align="right">Targeted</th>
    {}
    </table>
    """.format(data)


class HTMLoutput:
    """Context manager for simple HTML output of results
    """

    def __init__(self, html_fn, title):
        self.html_fn = os.path.abspath(html_fn)
        self.fig_dir = os.path.join(os.path.dirname(self.html_fn), 'figures')
        if not os.path.exists(self.fig_dir):
            os.makedirs(self.fig_dir)
        self.title = title
        self.fig_idx = 1

        self.summary = StringIO()
        self.details = StringIO()
        self.summary_metrics = []

        self.summary.write("""
          <h3>Summary Metrics</h3>
          <p>Summary metrics with targets.</p>
        """.format())

    def __enter__(self):
        self.outfile = open(self.html_fn, 'w')
        self.outfile.write("""
        <html><head>
        <title>scATAC report - {title}</title>
        <style>{css}</style>
        </head>
        <body>
        <h1>{title}</h1>

        <div class="tab">
          <button class="tablinks" onclick="openTab(event, 'Summary')" id="defaultOpen">Summary</button>
          <button class="tablinks" onclick="openTab(event, 'Details')">Details</button>
        </div>
        """.format(title=self.title, css=CSS))
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        if self.summary_metrics:
            metrics = sorted(self.summary_metrics, key=lambda m: m.category)
            self.summary.write(targeted_metric_list_to_html(metrics))

        self.outfile.write("""
        <div id="Summary" class="tabcontent">
            {summary}
        </div>

        <div id="Details" class="tabcontent">
            {details}
        </div>

        <script>{script}</script>
        <body></html>
        """.format(summary=self.summary.getvalue(),
                   details=self.details.getvalue(),
                   script=SCRIPT))
        self.outfile.close()

    def add_results(self, figures=None, title=None, metrics=None, text=None, targeted_metrics=None):
        """Adds a section of results, containing optionally figures, a table of metrics, and explanatory text"""
        if targeted_metrics is not None:
            self.summary_metrics.extend(targeted_metrics)

        if figures is None and metrics is None and text is None:
            # No detail section
            return

        if title is None:
            title = 'Section'

        self.details.write("""
        <button class="accordion">{}</button>
        <div class="panel">
        """.format(title))

        if text is not None:
            self.details.write('<p>{}</p>\n'.format(text))

        if metrics is not None:
            self.details.write(
                """<p class="metrics">{}</p>\n
                """.format(metric_list_to_html(metrics)))

        if figures is not None:
            for figure in figures:
                figout = os.path.abspath(os.path.join(self.fig_dir, 'figure{}.png'.format(self.fig_idx)))
                figure.savefig(figout, format='png')
                with open(figout, 'rb') as figdata:
                    encoded = base64.b64encode(figdata.read()).replace('\n', '')
                self.details.write('<p class="matplotlib"><img src="data:image/png;base64,{}"></p>\n'.format(encoded))
                self.fig_idx += 1

        self.details.write('</div>')
#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
# Test various utility functions

import tenkit.test as tk_test
from utils import is_url

class TestUtils(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_url_format_check(self):
        """Tests the code used to parse url and verify formats. Gives a slew of cases that are supposed to success as well as fail"""

        assert not is_url("www.google.com")
        assert is_url("http://www.google.com")
        assert is_url("https://www.google.com")
        assert is_url("http://www.google.com/info.txt")
        assert is_url("http://www.google.com/child/info.txt")

        assert not is_url("10.120.1.23")
        assert is_url("http://10.120.1.23")
        assert is_url("http://10.120.1.23/info.txt")
        assert is_url("http://10.120.1.23/child/info.txt")

        assert is_url("http://127.0.0.1:8080")
        assert is_url("http://127.0.0.1:8080/child/info.txt")

        assert is_url("http://port:8080")
        assert is_url("http://port:8080/child/info.txt")

        assert is_url("http://hello")
        assert not is_url("http://hello.")
        assert is_url("http://hello.i")
        assert is_url("http://hello.io")
        assert is_url("http://hello/child/info.txt")

        assert is_url("http://hel-lo")
        assert is_url("http://hel_lo")
        assert not is_url("http://hel lo")
        assert is_url("http://hello/")
        assert is_url("http://hello/.")
        assert is_url("http://hello/.txt")



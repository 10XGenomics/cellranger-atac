# Websummary


## Description of Files
* components.js - React components that gets compiled and minified
* template.html - template where the compiled components.js, vendored dependencies, data, and (websummary) javascript will get inlined
* summarize.py - simple script that inlines all the javascript into the template

## Expected use
call summarize.generate\_html\_summary() with a json-serializable object describing the data and the path
to an HTML template that uses the data. The example directory contains an example JSON bag and template.

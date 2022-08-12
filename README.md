# multi-methylprep

In order to make multi-processing idat files, replace below three files in a methylprep package.

- /home/user/.local/lib/python3.8/site-packages/methylprep/processing/pipeline.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/processing/postprocess.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/models/sigset.py

# postprocess.py
```
consolidate_values_for_sheets() -> consolidate_values() and consolidate_values_for_sheets()
```

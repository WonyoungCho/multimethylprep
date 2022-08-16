# multi-methylprep

In order to make multi-processing idat files, modify below three files in a methylprep package.

- /home/user/.local/lib/python3.8/site-packages/methylprep/processing/pipeline.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/processing/postprocess.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/models/sigset.py

# postprocess.py

`consolidate_values_for_sheets()` split into `consolidate_values()` and `consolidate_values_for_sheets()`.

```python
import parmap

def consolidate_values(batch_data_containers, postprocess_func_colname='beta_value', bit='float32', poobah=True, poobah_sig=0.05, exclude_rs=True, np=1):
    if np > 1:
        dfs = parmap.map(consolidate_values_for_sheets, batch_data_containers,
                         postprocess_func_colname=postprocess_func_colname, bit=bit,
                         poobah=poobah, poobah_sig=poobah_sig, exclude_rs=exclude_rs,
                         pm_pbar=True, pm_processes=np
                         )
        merged = pd.concat(dfs, axis=1)

    elif np == 1:
        merged = consolidate_values_for_sheet(batch_data_containers, postprocess_func_colname=postprocess_func_colname, bit=bit, poobah=poobah, exclude_rs=True)
        
    if bit != 'float32' and bit in ('float64','float16'):
        merged = merged.astype(bit)
        
    return merged


def consolidate_values_for_sheets(sample, postprocess_func_colname='beta_value', bit='float32', poobah=True, poobah_sig=0.05, exclude_rs=True):
    poobah_column = 'poobah_pval'
    quality_mask = 'quality_mask'

    sample_id = f"{sample.sample.sentrix_id}_{sample.sample.sentrix_position}"

    if poobah == True and poobah_column in sample._SampleDataContainer__data_frame.columns:
        # remove all failed probes by replacing with NaN before building DF.
        sample._SampleDataContainer__data_frame.loc[sample._SampleDataContainer__data_frame[poobah_column] >= poobah_sig, postprocess_func_colname] = np.nan
    elif poobah == True and poobah_column not in sample._SampleDataContainer__data_frame.columns:
        LOGGER.warning('DEBUG: missing poobah')

    if sample.quality_mask == True and quality_mask in sample._SampleDataContainer__data_frame.columns:
        # blank there probes where quality_mask == 0
        sample._SampleDataContainer__data_frame.loc[sample._SampleDataContainer__data_frame[quality_mask] == 0, postprocess_func_colname] = np.nan

    this_sample_values = sample._SampleDataContainer__data_frame[postprocess_func_colname]

    if exclude_rs: # dropping rows before exporting
        mask_snps = (sample._SampleDataContainer__data_frame.index.str.startswith('rs'))
        this_sample_values = this_sample_values.loc[ ~mask_snps ]

    df = pd.DataFrame(this_sample_values, columns=[postprocess_func_colname])
    df.rename(columns={postprocess_func_colname: sample_id}, inplace=True)

    return df
```

# sigset.py
This function reads `idat` file in the designated directory.

```python
import multiprocessing as mp
import parmap

def collect_idats(sample,zip_reader=None):
    green_filepath = sample.get_filepath('idat', Channel.GREEN)
    green_idat = IdatDataset(green_filepath, channel=Channel.GREEN)
    red_filepath = sample.get_filepath('idat', Channel.RED)
    red_idat = IdatDataset(red_filepath, channel=Channel.RED)
    return {'green_idat': green_idat, 'red_idat': red_idat, 'sample': sample}
    
def parse_sample_sheet_into_idat_datasets(sample_sheet, sample_name=None, from_s3=None, meta_only=False):
...
    if from_s3 and meta_only:
        parser = RawMetaDataset
        idat_datasets = [parser(sample) for sample in samples]
    elif from_s3 and not meta_only:
        #parser = RawDataset.from_sample_s3
        zip_reader = from_s3
        print('Reading IDATs ...')
        idat_datasets = parmap.map(collect_idats, samples, zip_reader, pm_pbar=True, pm_processes=mp.cpu_count())

    elif not from_s3 and not meta_only:
        #parser = RawDataset.from_sample
        print('Reading IDATs ...')
        idat_datasets = parmap.map(collect_idats, samples, pm_pbar=True, pm_processes=mp.cpu_count())
...
```

# pipeline.py
```python
from .postprocess import (
    ...
    consolidate_values,
    consolidate_values_for_sheet,
    ...
)

import parmap


```

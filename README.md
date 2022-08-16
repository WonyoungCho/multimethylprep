# multi-methylprep

In order to make multi-processing idat files, modify below three files in a methylprep package.

- /home/user/.local/lib/python3.8/site-packages/methylprep/processing/pipeline.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/processing/postprocess.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/models/sigset.py

# postprocess.py
```
consolidate_values_for_sheets() split into consolidate_values() and consolidate_values_for_sheets()
```

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


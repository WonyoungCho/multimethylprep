# multi-methylprep

```bash
$ pip install methylprep
```

In order to make multi-processing `methylprep (1.7.1)`, replace below four files in the methylprep package directory.

- /home/user/.local/lib/python3.8/site-packages/methylprep/cli.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/processing/pipeline.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/processing/postprocess.py
- /home/user/.local/lib/python3.8/site-packages/methylprep/models/sigset.py

# cli.py
```python
import multiprocessing as mp

def cli_process(cmd_args):
    ...
    parser.add_argument(
        '-th', '--threads',
        required=False,
        type=int,
        default=mp.cpu_count(),
        help='Number of threads to run jobs (default : '+str(mp.cpu_count())+' ).'
    )
    ...
    
    run_pipeline(
        ...
        np=args.threads    
        ...
```

# sigset.py
- This function reads `idat` file in the designated directory.

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

# postprocess.py
- Insert `consolidate_values()` and `consolidate_values_for_sheets()`.
- Change parquet suffix `parquet` to `par`.

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

# pipeline.py

- Replace `consolidate_values_for_sheet()` to `consolidate_values(..., np=np)` in `if {beta_values, m_values, ...}`.
- Take `_prepare_save_out_file()` out and add more arguments, `_prepare_save_out_file(data_dir, df, file_stem, batch_size, batch_num, file_format, uint16=False)`.
- Change parquet suffix `parquet` to `par`.

```python
from .postprocess import (
    ...
    consolidate_values,
    consolidate_values_for_sheet,
    ...
)

import parmap

def processing_idats(idat_dataset_pair, manifest, save_uncorrected, bit, do_infer_channel_switch,
                     sesame, quality_mask, do_noob, poobah, poobah_decimals, poobah_sig,
                     do_nonlinear_dye_bias, kwargs,
                     export, save_control, low_memory):
    data_container = SampleDataContainer(
        idat_dataset_pair,
        manifest=manifest,
        retain_uncorrected_probe_intensities=save_uncorrected,
        bit=bit,
        switch_probes=(do_infer_channel_switch or sesame), # this applies all sesame-specific options
        quality_mask= (quality_mask or sesame or False), # this applies all sesame-specific options (beta / noob offsets too)
        do_noob=(do_noob if do_noob != None else True), # None becomes True, but make_pipeline can override with False
        pval=poobah, #defaults to False as of v1.4.0
        poobah_decimals=poobah_decimals,
        poobah_sig=poobah_sig,
        do_nonlinear_dye_bias=do_nonlinear_dye_bias, # start of run_pipeline sets this to True, False, or None
        debug=kwargs.get('debug',False),
        sesame=sesame,
    )

    data_container.process_all()

    output_path = None;
    noob_error = None; raw_error = None;
    sample_id = None; control_df = None

    if export: # as CSV
        output_path = data_container.sample.get_export_filepath()
        data_container.export(output_path)

        # this tidies-up the tqdm by moving errors to end of batch warning.
        if data_container.noob_processing_missing_probe_errors != []:
            noob_error = data_container.noob_processing_missing_probe_errors

        if data_container.raw_processing_missing_probe_errors != []:
            raw_error = data_container.raw_processing_missing_probe_errors

    if save_control: # Process and consolidate now. Keep in memory. These files are small.
        sample_id = f"{data_container.sample.sentrix_id}_{data_container.sample.sentrix_position}"
        control_df = one_sample_control_snp(data_container)

    # now I can drop all the unneeded stuff from each SampleDataContainer (400MB per sample becomes 92MB)
    # these are stored in SampleDataContainer.__data_frame for processing.
    if low_memory is True:
        # use data_frame values instead of these class objects, because they're not in sesame SigSets.
        del data_container.man
        del data_container.snp_man
        del data_container.ctl_man
        del data_container.green_idat
        del data_container.red_idat
        del data_container.data_channel
        del data_container.methylated
        del data_container.unmethylated
        del data_container.oobG
        del data_container.oobR
        del data_container.ibG
        del data_container.ibR

    return data_container, sample_id, control_df, noob_error, raw_error, output_path

    
def _prepare_save_out_file(data_dir, df, file_stem, batch_size, batch_num, file_format, uint16=False):
    out_name = f"{file_stem}_{batch_num}" if batch_size else file_stem
    if uint16 and file_format != 'parquet':
        df = df.astype('float32') if df.isna().sum().sum() > 0 else df.astype('uint16')
    else:
        df = df.astype('float32')
    if df.shape[1] > df.shape[0]:
        df = df.transpose() # put probes as columns for faster loading.
    # sort sample names
    df = df.sort_index().reindex(sorted(df.columns), axis=1)
    if file_format == 'parquet':
        # put probes in rows; format is optimized for same-type storage so it won't really matter
        df.to_parquet(Path(data_dir,f"{out_name}.par"))
    else:
        df.to_pickle(Path(data_dir, f"{out_name}.pkl"))
    LOGGER.info(f"saved {out_name}")
    
    
def run_pipeline(np=1,):
    ...
    for batch_num, batch in enumerate(batches, 1):
    ...
        batch_data_containers = []
        export_paths = set() # inform CLI user where to look

        print('Processing samples ...')
        data_containers = parmap.map(
            processing_idats, idat_datasets, manifest, save_uncorrected, bit,
            do_infer_channel_switch, sesame, quality_mask, do_noob, poobah,
            poobah_decimals, poobah_sig, do_nonlinear_dye_bias, kwargs,
            export, save_control, low_memory,
            pm_pbar=True, pm_processes=np
        )

        del idat_datasets

        for i in range(len(data_containers)):
            batch_data_containers.append(data_containers[i][0])
            if data_containers[i][1] != None: control_snps[data_containers[i][1]] = data_containers[i][2]
            if data_containers[i][3] != None: missing_probe_errors['noob'].extend(data_containers[i][3])
            if data_containers[i][4] != None: missing_probe_errors['raw'].extend(data_containers[i+4])
            if data_containers[i][5] != None: export_paths.add(data_containers[i][5])

        del data_containers

        if kwargs.get('debug'): LOGGER.info('[finished SampleDataContainer processing]')

        funcDict = {'beta_value':'beta','m_value':'m','noob_meth':'noob_meth','noob_unmeth':'noob_unmeth',
                    'meth':'meth', 'unmeth':'unmeth', 'poobah_pval':'poobah'}

        if betas:
            df = consolidate_values(batch_data_containers, postprocess_func_colname='beta_value', bit=bit, poobah=poobah, exclude_rs=True, np=np)
            _prepare_save_out_file(data_dir, df, 'beta_values', batch_size, batch_num, file_format)

        if m_value:
            df = consolidate_values(batch_data_containers, postprocess_func_colname='m_value', bit=bit, poobah=poobah, exclude_rs=True, np=np)
            _prepare_save_out_file(data_dir, df, 'm_values', batch_size, batch_num, file_format)

        if (do_save_noob is not False) or betas or m_value:
            df = consolidate_values(batch_data_containers, postprocess_func_colname='noob_meth', bit=bit, poobah=poobah, exclude_rs=True, np=np)
            _prepare_save_out_file(data_dir, df, 'noob_meth_values', batch_size, batch_num, file_format, uint16=True)

            # TWO PARTS
            df = consolidate_values(batch_data_containers, postprocess_func_colname='noob_unmeth', bit=bit, poobah=poobah, exclude_rs=True, np=np)
            _prepare_save_out_file(data_dir, df, 'noob_unmeth_values', batch_size, batch_num, file_format, uint16=True)

        #if (betas or m_value) and save_uncorrected:
        if save_uncorrected:
            df = consolidate_values(batch_data_containers, postprocess_func_colname='meth', bit=bit, poobah=False, exclude_rs=True, np=np)
            _prepare_save_out_file(data_dir, df, 'meth_values', batch_size, batch_num, file_format, uint16=True)

            # TWO PARTS
            df = consolidate_values(batch_data_containers, postprocess_func_colname='unmeth', bit=bit, poobah=False, exclude_rs=True, np=np)
            _prepare_save_out_file(data_dir, df, 'unmeth_values', batch_size, batch_num, file_format, uint16=True)

        if export_poobah:
            if all(['poobah_pval' in e._SampleDataContainer__data_frame.columns for e in batch_data_containers]):
                # this option will save a pickled dataframe of the pvalues for all samples, with sample_ids in the column headings and probe names in index.
                # this sets poobah to false in kwargs, otherwise some pvalues would be NaN I think.
                df = consolidate_values(batch_data_containers, postprocess_func_colname='poobah_pval', bit=bit, poobah=False, poobah_sig=poobah_sig, exclude_rs=True, np=np)
                _prepare_save_out_file(data_dir, df, 'poobah_values', batch_size, batch_num, file_format)
                
            if all(['pNegECDF_pval' in e._SampleDataContainer__data_frame.columns for e in batch_data_containers]):
                # this option will save negative control based pvalues for all samples, with
                # sample_ids in the column headings and probe names in index.
                df = consolidate_values(batch_data_containers, postprocess_func_colname='pNegECDF_pval', bit=bit, poobah=False, poobah_sig=poobah_sig, exclude_rs=True, np=np)
                _prepare_save_out_file(data_dir, df, 'pNegECDF_values', batch_size, batch_num, file_format)
    ...
```

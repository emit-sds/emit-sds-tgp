<div align="center">

# Trace Gas Plume Detection and Quantification

</div>


This repository holds the control mechanism to interface with manual or ML based delineations, and adds metadata, filtering, quantifications, and delivery support. There are several core breakouts this repository:

- Documentation, including the [Algorithm Theoretical Bases Document](docs/EMIT_L2B_TRACE_GAS_ATBD.md) and [User Guide](docs/EMIT_L2B_TRACE_GAS_User_Guide.md)
- Annotation support, which is mainly for use by the EMIT SDS team
- [Quantification](src/quantification/quantification_readme.md), much of which can be run independently of the rest of the workflow, which may be of interest outside of operational workflow

The general control flow used in operations is scrape_refine_upload, executed as:

```
python src/scrape_refine_upload.py {token} {layer} {output_dir} 002 --loglevel INFO --database_config /store/emit/ops/repos/emit-main/emit_main/config/ops_sds_config.json --track_coverage_file /store/brodrick/emit/emit-visuals/track_coverage_pub.json
```



# Background

The configuration information for the workflow is contained in "flat" key-value files
like the following:

```
key1="value1"
key2="value2"
...
keyN="value3"
```

However, Cromwell/WDL takes a JSON format file like the following:

```
{
    "key1": "value1",
    "key2": "value2",
    ...
    "keyN": "valueN"
}

```

This parser takes one or more files in the "flat" key-value format, and converts them into a JSON 
  formatted file.
  
# Details

## Index files

Many index files are implicitly present in the same location as the primary file, such as

```
# (Main reference file)
ref.fasta

# (Index files for bwa)
ref.fasta.ann
ref.fasta.bwt
...
etc.
```

Often, these index files are not passed explicitly to the tools that use the main file, 
  but they are still required for the tool to run properly.

To ensure we could validate that all inputs exist, we decided that implicit files should
  be made explicit in the JSON file.
  
Therefore, the parser treats variables like "Ref" and "DBSNP" as special keys, and automatically
  adds their index files to the output JSON file.
  
  For example, if the flat file is:
  
```
DBSNP=dbsnp.vcf
```
The JSON file will have the following:
```
{
    "DBSNP": "dbsnp.vcf",
    "DBSNPIdx: "dbsnp.vcf.idx"
}
```

## Basic checks

While type validation occurs on the JSON formatted file that is produced by the parser, some
  basic checks are done during parsing.
  
Basic checks:

    1. Is every value surrounded by quote marks?
    2. Are any values empty strings?
    3. Are any keys present multiple times?
    4. Are any values completely whitespace characters surrounded by quote marks?
    5. Do any values have special characters? (Any of these: !#$%&()*;<>?@[]^`{|}~ )
    
For values found in the 'OPTIONAL_KEYS' list in 'src/config/util/special_keys.py', these values are
  allowed to be empty. If they are, these basic checks are not done 
## adding subworkflows (already done)
```
git submodule add git@github.com:oicr-gsi/mutect2.git subworkflows/mutect2
git submodule add git@github.com:oicr-gsi/variantEffectPredictor.git subworkflows/variantEffectPredictor
git submodule add git@github.com:oicr-gsi/HSMetrics.git subworkflows/HSMetrics
```

## update subworkflow repos
```
git submodule update --init --recursive
```

Or git checkout the desired subworkflow submodule commit/tag version.

## update subworkflow "pull" imports WDL

This step "pulls up" all subworkflow task level parameters into a workflow level parameters so that they are addressable from the primary workflow.

```
module load gsi-wdl-tools

generate-subworkflow-import --input-wdl subworkflows/HSMetrics/hsMetrics.wdl --pull-all --output-wdl-path imports/pull_hsMetrics.wdl
generate-subworkflow-import --input-wdl subworkflows/mutect2/mutect2.wdl --pull-all --output-wdl-path imports/pull_mutect2.wdl
generate-subworkflow-import --input-wdl subworkflows/variantEffectPredictor/variantEffectPredictor.wdl --pull-all --output-wdl-path imports/pull_variantEffectPredictor.wdl
```

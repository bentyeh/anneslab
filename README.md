# STF-1081 Drug Toxicity Target Identification

## Objective

Identify off-target kinases inhibited by STF-1081 that are likely responsible for its toxicity.

## Background

OTSSP167 is a potent MELK inhibitor (IC<sub>50</sub> = 0.41 nM) developed as a potential therapeutic for a wide range of cancers. [[Chung]](#ref) Of interest to our lab, it also potently inhibits DYRK1A, a potential target for pancreatic beta-cell replication drugs. [[Wang]](#ref) However, *in vivo* mice studies (citation needed) suggest that OTSSP167 is too toxic to be repurposed directly as a beta-cell replication drug. Our data supports the hypothesis that off-target promiscuity is in part responsible for this toxicity: 138 out of 403 non-mutant kinases in the KINOMEscan panel were inhibited to below 10% activity by 100 nM OTSSP167.

Rather than perform a costly systematic or high-throughput knockdown/knockout screen of all inhibited targets of OTSSP167, we used data from our own lab and published literature to identify likely candidates responsible for its toxicity.

## Datasets

Datasets were chosen from the literature based on inclusion of both OTSSP167 and one or more less-toxic compounds against a broad kinase panel.

| Name    | Kinase panel                                     | # targets tested | Metric                  | Compounds tested (concentration)             |
| ------- |:------------------------------------------------ |:---------------- |:----------------------- |:-------------------------------------------- |
| Annes   | KINOMEscan scanMAX                               | 468 \*\*         | % control \*            | OTSSP167 (100 nM), STF-1285 (100 nM, 500 nM) |
| Klaeger | Kinobeads                                        | 520 \*\*\*       | K<sub>d,app</sub> \*    | OTSSP167, CC-401, many others (3 nM - 30 uM) |
| Huang   | International Centre for Kinase Profiling (ICKP) | 140              | % activity remaining \* | OTSSP167, HTH-01-091 (1 uM)                  |

\* Lower value indicates more effective inhibition.
\*\* Includes kinases with multiple variants (different phosphorylation states or mutations)
\*\*\* Includes 242 protein and lipid kinases, 88 non-kinase direct binders ("nucleotide binders, helicases, ATPases and GTPases, FAD (e.g., NQO2) and heme (e.g., FECH) containing proteins"), and 190 non-kinase non-direct binders ("interaction partners/adaptor proteins of the kinases"). See the "Target selection criteria" section within the *Supplementary Text* chapter under Supplementary Materials for Klaeger et al.

## Methods <a name="methods"></a>

### Boolean filtering
Find the intersection across all datasets of targets that are inhibited by OTSSP167 but not by a less toxic compound, such as CC-401 [Celgene], HTH-01-091 [[Huang]](#ref), or STF-1285 [Annes].

### Rank ordering
Within each dataset, rank each target by the difference in inhibition between OTSSP167 and a less toxic drug; normalize by the number of targets in the dataset. For each target, sum the rank across all datasets. Lower rank values indicate more differentially inhibited by OTSSP167 than the less toxic compounds.

## Files

data/: raw data
* Annes100.csv, Annes500.csv
  * Concatenated raw data from DiscoverX KINOMEscan at 100 nM OTSSP167 and 100 or 500 nM STF-1285.
* Huang.csv
  * Downloaded as a [Microsoft Word file](https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMjY2OTMvZWxpZmUtMjY2OTMtZmlnMS1kYXRhMS12MS5kb2N4/elife-26693-fig1-data1-v1.docx?_hash=UbxyMBf8pkVSYVLfuzWLRfqv35PCZ3rDoDex%2B6QmdWc%3D) from [Figure 1---source data 1](https://doi.org/10.7554/eLife.26693.004). Copied only the "Kinase" and "% AR" columns to a CSV spreadsheet. Renamed columns as "Name", "OTSSP167", and "HTH01091" so that column names are syntactically valid in R.
* Klaeger.csv
  * Extracted "Kinobeads Drugmatrix detailed" sheet within aan4368_Table_S3.xlsx downloaded from [Supplementary Materials](http://science.sciencemag.org/highwire/filestream/702936/field_highwire_adjunct_files/1/aan4368_Tables1-11.zip) of [Huang et al.](#ref). Kept only "Name", "CC-401" and "OTS-167" columns. Renamed "CC-401" and "OTS-167" to "CC401" and "OTSSP167" to so that column names are syntactically valid in R.
  * Note: If using Microsoft Excel to view this file, low-confidence values in parentheses may appear as negative values, even though the parentheses are properly displayed in a normal text editor and are kept when reading the CSV file using standard packages such as pandas (Python) or readxl (R). See this [StackOverflow](https://stackoverflow.com/questions/29648572/excel-values-in-parentheses-become-negative) post for a discussion about how/why Excel handles parenthesized numbers.
    * Definition of high/low-confidence values according to [Klaeger, et al.](#ref): "A protein was considered a high-confidence target if the binding curve showed a sigmoidal shape with a dose-dependent decrease in binding to the Kinobeads."

data_aux/: processed data files produced by scripts/process_data.ipynb
* Annes100.csv, Annes500.csv
  * For kinases with multiple variants (different phosphorylation states, mutations), only keep the mean value.
  * Convert target names to HGNC official gene symbols.
* Huang.csv: Convert kinase names to HGNC official gene symbols.
  * Convert target names to HGNC official gene symbols.
* Klaeger.csv: 
  * Convert target names to HGNC official gene symbols.

scripts/:
* process_data.ipynb
  * Process data as described for each file in data_aux/ above
* selectTargets.ipynb
  * Select potential targets as described in [Methods](#methods) section

results/: results files produced by scripts/selectTargets.ipynb
* intersect.txt
  * List of targets identified using the Boolean filtering method.
* rank.tsv
  * Targets ranked according to the Rank ordering method.

### Dependencies

Python 3
- Biopython
- pandas
- numpy
- (optional) Jupyter - to run the code in a Jupyter notebook
- (optional) IPython - to run the code in a Jupyter notebook

## Results

Initial *in vitro* validation experiments suggest that Aurora kinase B (AURKB) is in part responsible for the toxicity of OTSSP167.

## References <a name="ref"></a>
1. [Klaeger, S. et al. The target landscape of clinical kinase drugs. *Science* 358, eaan4368 (2017).](http://science.sciencemag.org/content/358/6367/eaan4368)
2. [Huang, H.-T. et al. MELK is not necessary for the proliferation of basal-like breast cancer cells. *Elife* 6, e26693 (2017).](https://elifesciences.org/articles/26693)
3. [Chung, S. et al. Development of an orally-administrative MELK-targeting inhibitor that suppresses the growth of various types of human cancer. *Oncotarget* 3, 1629-40 (2012)](https://doi.org/10.18632/oncotarget.790)
4. [Wang, P. et al. A high-throughput chemical screen reveals that harmine-mediated inhibition of DYRK1A increases human pancreatic beta cell replication. *Nature Medicine* 21, 383-388 (2015).](https://www.nature.com/articles/nm.3820)
5. [DiscoverX KINOMEscan(R)](https://www.discoverx.com/services/drug-discovery-development-services/kinase-profiling/kinomescan)
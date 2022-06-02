# Expression explorer

## Demo

You can fire-up a demo of this application on binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/basf/expression_explorer/master?urlpath=shiny/app/)

## Deploy with docker (local or server)

Requires: docker-compose

Run: docker-compose up --build -d
The expression explorer will then be available at:
http://localhost:8005/expression_explorer/

# Adding a dataset

Each dataset resides in its dedicated folder, under "data". The name of the folder will be displayed in the starting table, together with additional information about the dataset. Adding actual data is as easy as creating one of the standard and simple tab separated files that EE expects. Some of these are required for the processing of the dataset, others are optional.

An example dataset can be found under "example". Once you deploy the app to your favourite server, you can simply delete this folder.

## Gene Counts
Raw, non-normalized read counts should be provided in the "counts.tab" file. This file should have unique gene identifiers as the first column, and per-sample counts for any addititional column. The header of this file must contain the repective sample names.
For example:
| gene | sample1 | sample2 | sample3 |
| ------ | ------ | ------ | ------ |
| gene1 | 1 | 2 | 2 |
| gene2 | 4 | 4 | 2 |
| gene3 | 2 | 6 | 5 |

## Gene annotations
For each gene identifier, additional information/annotations can be provided in the "ann.tab" file. File header is as follows:
| gene | annotation | synonym |
| ------ | ------ | ------ |

The synonym column is the one used when plotting information about the genes, so we recommend using short, readable names in this column. The annotation should be a long, descriptive one. Note that this annotation can, by default, not be edited by the user of the app. It is however possible to alow the user to edit this information.

## Medadata
Sample metadata is very important for making sense of any dataset. EE requires four columns, but allows to user to add as much metadata as they please. This information should be present in the "meta.tab" file. The first columns and the mandatory header are:
| sample_name | Sample | Timepoint | Replicate |
| ------ | ------ | ------ | ------ |

The sample_name must match the sample names from the "counts.tab" file.

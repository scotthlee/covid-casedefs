# Household transmission case definition study

## Overview
Describe the purpose of your project. Add additional sections as necessary to help collaborators and potential collaborators understand and use your project.

## The code
We used both R and Python to do the analysis for this study. There's a good deal of overlap between the two sets of code, but in general, the R code is geared toward generating our graphics, while the Python code is geared toward generating the statistics we report in the tables. Here's quick manifest of the scripts:

  ### Python
  1. (tools.py)[https://github.com/scotthlee/hh-transmission/blob/master/python/tools.py]: support functions used for the analysis
  2. (multi.py)[https://github.com/scotthlee/hh-transmission/blob/master/python/multi.py]: multiprocessing-enabled versions of functions from `tools.py`
  3. (combo_search.py)[https://github.com/scotthlee/hh-transmission/blob/master/python/combo_search.py]: runs the combinatorial symptom search
  4. (rf.py)[https://github.com/scotthlee/hh-transmission/blob/master/python/rf.py]: trains a random forest on the data
  5. (primary_analysis.py)[https://github.com/scotthlee/hh-transmission/blob/master/python/primary_analysis.py]: produces the statistics and tables in the manuscript


For more info, check out the respective READMEs. 

## The data
In the data, we have information about symptoms for each of the study participants, and we also have their SARS-CoV-2 PCR and ELISA test results. In our primary analysis, we frame the problem as one of binary classification, i.e., by predicting PCR status from different combinations of the symptoms. Here's a quick rundown of the variables:

  1. `study_id`: participant identifier
  2. `hh_id`: household identifier
  3. `age_adult`: whether age is over (1) or under (0) 18 years
  4. `wheeze` to `tastesmell_combo`: the symptoms
  5. `ili`: influenza-like illness
  6. `cdc`: CDC's COVID symptom list
  7. `ari`: the WHO RSV ARI case definition
  8. `cste`: the CSTE COVID case definition
  9. `cli`: COVID-like illness
  10. `vaccine`s: CDC COVID vaccine group proposed trial endpoints
  11. `sero_pos`: whether ELISA detected SARS-CoV-2 antibodies
  12. `sero_conv`: whether the contact seroconverted during the 2-week observation period
  13. `pcr_pos`: whether RT-PCR detected SARS-CoV-2 infection
  14. `any_pos`: whether ELISA or RT-PCR was positive

## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/privacy.html](http://www.cdc.gov/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page are subject to the [Presidential Records Act](http://www.archives.gov/about/laws/presidential-records.html)
and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).

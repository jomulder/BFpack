# BFpack 1.1.0

* Date: 2022-10-31
* Extension of BF.rma.uni for also allowing a fixed effects meta-analytic model.
* Including an additional argument BF.type for different Bayes factor tests for the same model. Currently only BF.type = 1 and 2 are support for a model class 'lm' or 'mlm' or a named vector where BF.type = 1 (implies no prior adjustment to the null value) and BF.type = 2 implies a prior adjustment to the null value. The default choice is BF.type = 2.
* Bugs fixes


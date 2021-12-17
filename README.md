# ACD_methods
Code and documentation to perform automated change detection methods using Google Earth Engine.

## Active repository

This is an active repository currently has Earth Engine JavaScript files for:

- Performing automated change detection following the general approach used to update the [National Land Cover Dataset](http://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=1720&context=usgsstaffpub&sei-redir=1&referer=https%3A%2F%2Fscholar.google.com%2Fscholar%3Fhl%3Den%26as_sdt%3D0%252C9%26q%3Dnlcd%2Bland%2Bcover%2Bchange%26btnG%3D#search=%22nlcd%20land%20cover%20change%22)
- Performing automated change detection using iteratively reweighted multivariate alteration detection (IW-MAD), as described by [Nielsen, A.A. 2007](https://pdfs.semanticscholar.org/ec15/ae9a4b26f07d0ba9e0c952daf58cae8b243e.pdf)
- Masking clouds from Sentinel 2 images.
- Calculating and appending shape metrics to polygon attributes.

## Usage

Currently, script provided in the EE folder will need to be loaded into and run in Google Earth Engine code editor.

To see an example of an analysis using these methods, see our [case study](http://cci-dev.org/analysis/LPC_delisting) detecting energy development in the range of the Lesser prairie chicken.

For detailed explanation of methodlology, see our [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/611459v1)

The Algorithm Validation directory contains R code and files to perform LDA and AUC analyses to assess the sensitivity and specificity of algorithm outputs based on training and testing data.

## Contributing

Interested in applying or helping improve change detection methods? [Get in touch with us](mailto:esa@defenders.org).

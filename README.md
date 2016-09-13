# correlation_analysis

-----------------------------------------------------------------------
Metabolites Correlation Analysis Tool
-----------------------------------------------------------------------

This package contains the metabolite_correlation_analysis tool version="20141118".


--------------------------------------------------------------------
Tests on R bin and Packages
--------------------------------------------------------------------

$ R --version
R version 2.15.1 (2012-06-22) -- "Roasted Marshmallows"
Platform: x86_64-unknown-linux-gnu (64-bit)

R packages needed for this tool:

library(batch) #necessary for parseCommandArgs function
library(reshape) #necessary for using melt function
library(MASS) # necessary for using the write.matrix()


--------------------------------------------------------------------
Instructions for integration of the this tool into the workflow-system
Galaxy (http://getgalaxy.org)
--------------------------------------------------------------------

For installing the tools of the XCMS Suite into your Galaxy installation, please do the following:

	- First of all, you need to have an account on the ABIMS server. Follow the link to abims website, where it is explained how to get an account: 
http://abims.sb-roscoff.fr/faq

	- For a manual installation, copy the image from the  static/images/ into the  "galaxy-dist/static/images/"



Last but not least, restart Galaxy.


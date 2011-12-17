\name{pcaGoPromoter-package}
\Rdversion{1.1}
\alias{pcaGoPromoter-package}
\alias{pcaGoPromoter}
\docType{package}
\title{
Tools for analyzing data from DNA micro arrays.
}
\description{
Welcome to the pcaGoPromoter package. The purpose of this package is to ease the analyses of data from genome wide gene expression analysis. Multiple platforms can be used (Affymetrix GeneChip, Illumina expression beadchips, qRT-PCR and more). Probe identifiers should be either Affymetrix Probe Set ID, Gene Symbol or Entrez gene ID.

The overall analysis strategy is to conduct principal component analysis followed by the interpretation of the principal component dimensions by functional annotation for biological function and prediction of regulatory transcription factor networks. 

You can learn what objects this package supports with the following command:

     'ls("package:pcaGoPromoter")' 

Each of these objects has their own manual page detailing where relevant data was obtained along with some examples of how to use it.
}
\author{
Morten Hansen \email{mhansen@sund.ku.dk} and Jorgen Olsen \email{jolsen@sund.ku.dk}

Maintainer: Morten Hansen \email{mhansen@sund.ku.dk}
}
\keyword{ package }
\examples{
ls("package:pcaGoPromoter")
}

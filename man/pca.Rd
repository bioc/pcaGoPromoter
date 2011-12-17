\name{pca}
\Rdversion{1.1}
\alias{pca}
\alias{print.pca}
\alias{plot.pca}
\alias{getRankedProbeIds}
\alias{getRankedProbeIds.pca}
\alias{pcaInfoPlot}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pca - Principal Component Analysis (PCA)
}
\description{
pca - This is a wrapper function around \code{prcomp} from stats package. It
will prepare the DNA micro array data, run the principal component analysis and
return an object of type 'pca'

plot.pca - Makes a score plot with legend of object type 'pca'.

getProbeIds.pca - Get probe ids from PCA loadings.

pcaInfoPlot - Makes a informative score plot with optional overrepresented GO
terms and predicted transcription factor binding sites annotated along the axis.
}
\usage{
pca(exprsData, printDropped = TRUE, scale=TRUE, center=TRUE)

\method{plot}{pca}(x, groups, PCs = c(1,2), printNames = TRUE,
     symbolColors = TRUE, plotCI = TRUE, GOtreeObjs, primoObjs,
     main, ...)

\method{getRankedProbeIds}{pca}(x, pc = 1, decreasing = TRUE)

pcaInfoPlot(exprsData, inputType = "hgu133plus2", org = "Hs", groups,
            noProbes = 1365, GOtermsAnnotation = TRUE,
            primoAnnotation = TRUE)
}
\arguments{
  \item{exprsData}{
a table with expression data. Variables in rows and observations in columns.
Row names should be probe identifiers (Affymetrix probe set ID, Gene Symbol or
Entrez IDs), column names should be sample identifiers.
}
  \item{printDropped}{
a logical value determination whether or not dropped probes should be printed.
Probes are dropped if all values are equal.
}
  \item{scale}{
a logical value determination whether or not input values should be scaled.
Default is TRUE.
}
  \item{center}{
a logical value determination whether or not input values should be centered.
Default is TRUE.
}
  \item{x}{
Object of type 'pca'. Created by function \code{pca}.
}
  \item{groups}{
a factor containing group information.
}

  \item{PCs}{
an integer vector indication which principal components to use in the score
plot. Exactly two principal components should be given.
}
  \item{printNames}{
a logical value indicating whether observation names should be printed.
Default is TRUE.
}
  \item{symbolColors}{
a logical value indicating whether the symbols should be plotted with colors.
Default is TRUE.
}
  \item{plotCI}{
a logical value indicating whether the confidence intervals for each group
should be plotted. Default is TRUE.
}

  \item{GOtreeObjs}{
a list of GOtree objects which is used for annotation of the axis in the PCA
plot - used by function pcaInfoPlot().
}
  \item{primoObjs}{
a list of primo objects which is used for annotation of the axis in the PCA
plot - used by function pcaInfoPlot().
}
  \item{main}{
a character vector with the title of the plot. If no title is given, a default
one will be provided. Set to NULL for no title.
}
  \item{...}{
other parameters to be passed through to plotting functions.
}
  \item{pc}{
a number indication which principal component to extract the probe identifiers
from.
}
  \item{decreasing}{
a logical value indication whether the probe identifiers from PCA loadings 
should be sorted in decreasing of ascending order
( \code{decreasing == FALSE} ).
}
  \item{inputType}{
a character vector with the input type. See '?primo' for details. Default is
"hgu133plus2".
}
  \item{org}{
a character vector with the organism. See '?primo' for details. Default is "Hs".
}
  \item{noProbes}{
a numeric value indication the number of probes to use for GOtree and primo.
}
  \item{GOtermsAnnotation}{
a logical value to select GO terms annotation along the axis. The default is
TRUE.
}
  \item{primoAnnotation}{
a logical value to select predicted transcription factor binding site
annotation along the axis. The default is TRUE.
}
}
\details{
pca uses \code{prcomp} to do the principal component analysis. The input data
(x) is scaled and centered, so constant variables (sd = 0) will be removed to
avoid divison by zero.

plot.pca makes a PCA score plot. The input data is an object of type 'pca'. The
score plot is formatted and the observations are colored according to the class.
A legend describing the class is placed below the plot. The proportion of
variance is plotted along the axis.

getRankedProbeIds.pca extract the probe identifiers from PCA loadings from an
object of type \code{pca}. The purpose of this function is to ease the
extraction of probe identifiers from an object of type \code{pca}. It generally
just gets the rownames from the loadings and sort it according to \code{pc} and
\code{decreasing}.
}
\value{
\code{pca} returns the an object of type 'pca'.

\code{getRankedProbeIds} return a vector with the probe ids sorted from higher
to lower if \code{decreasing = TRUE} or from lower to higher if
\code{decreasing = FALSE}.
}
\author{
Morten Hansen \email{mhansen@sund.ku.dk} and Jorgen Olsen
\email{jolsen@sund.ku.dk}
}
\seealso{
\code{prcomp}.
}
\examples{
library(serumStimulation)
data(serumStimulation)
groups <- as.factor(c(rep("control", 5), rep("serumInhib", 5),
                    rep("serumOnly", 3)))
pcaOutput <- pca(serumStimulation)
pcaOutput
plot(pcaOutput, groups=groups)
posLoadings <- getRankedProbeIds(pcaOutput)
}
\keyword{ methods }

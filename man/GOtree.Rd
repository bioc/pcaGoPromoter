\name{GOtree}
\Rdversion{1.1}
\alias{GOtree}
\alias{print.GOtree}
\alias{plot.GOtree}
\alias{GOtreeHits}
\alias{GOtreeWithLeaveOut}
\title{
GOtree, plot.GOtree and GOtreeWithLeaveOneOut
}
\description{
GOtree finds significantly overrepresented Gene Ontology terms in a list of
probes (only Biological processes) and return an object of type 'GOtree'.

print.GOtree list significant GO terms from an object of type 'GOtree'.

plot.GOtree creates a visual representation of the GO connection from an object
of type 'GOtree'.

GOtreeHits return the genes/probes for a specific GO term.

GOtreeWithLeaveOut returns the same as \code{GOtree}, but run through the
samples multiple times with 'Leave one out' cross-validation.
}
\usage{
GOtree(input, inputType = "hgu133plus2", org = "Hs",
       statisticalTest = "binom", binomAlpha = NA,
       p.adjust.method = "fdr")

\method{print}{GOtree}(x, ...)

\method{plot}{GOtree}(x, boxes = 25, legendPosition = "topright",
     main = "Gene Ontology tree, biological processes", ...)

GOtreeHits(input, inputType = "hgu133plus2", org = "Hs",
           GOid, returnGeneSymbols = TRUE)

GOtreeWithLeaveOut(exprsData, inputType = "hgu133plus2", org = "Hs",
                   pc = 1, decreasing = TRUE, noProbes = 1000,
                   leaveOut = 1, runs = NCOL(exprsData))
}
\arguments{
  \item{input}{
a character vector of Affymetrix probe ids, gene symbols or Entres gene IDs.
}
  \item{inputType}{
a character vector description the input type. Must be Affymetrix chip type,
"geneSymbol" or "entrezID". The following Affymetrix chip type are supported:
hgu133plus2, mouse4302, rat2302, hugene10st and mogene10st. Default is
Affymetrix chip type "hgu133plus2".
}
  \item{org}{
a character vector with the organism. Can be "Hs", "Mm" or "Rn". Only needed if
\code{inputType} is "geneSymbol" or "entrezID". See details. Default is "Hs".
}
  \item{statisticalTest}{
a character vector with the statistical method to be used. Can be "binom" or
"fisher". Default is "binom".
}
  \item{binomAlpha}{
a value with the pvalue for use in self contained test.
}
  \item{p.adjust.method}{
the method for adjust p-values due to multiple testing. This will come in a
separate column.
}
  \item{x}{
an object of type 'GOtree'. 
}
  \item{boxes}{
an integer indication the amount of boxes (terms) in the plot.
}
  \item{legendPosition}{
a vector description the position of the legend. See ?xy.coords for
possibilities. Set to NULL for no legend. Default is "topright".
}
  \item{main}{
a title for the GO tree plot
}
  \item{...}{
other parameters to be passed through to plotting functions.
}
  \item{GOid}{
a vector with the GO term of interest.
}
  \item{returnGeneSymbols}{
a logical indication whether gene symbols or probe ids should be returned.
Default is gene symbols.
}
  \item{exprsData}{
A table with expression data. Row names should be probe identifiers (Affymetrix
Probe set ID, Gene Symbols or Entrez gene ID). Column names should be sample
identifiers.
}
  \item{pc}{
a number indication which principal component to extract the probe list based
on the loading values from the pca.
}
  \item{decreasing}{
a logical value indication whether the loadings should be sorted in decreasing
of ascending order ( \code{decreasing == FALSE} ). Decreasing order yields
information about the positive direction and ascending order about the negative
direction of the particular principal component.
}
  \item{noProbes}{
a number indicating the number of probes included in the calculations
}
  \item{leaveOut}{
a number indication what percentage to leave out in the cross-validation. If
set to 1, each observation would be left out once and runs is set equal to
number of observations. Deafault is 1.
}
  \item{runs}{
a number indicating how many times to run with leave out. If leaveOut = 1, runs
is overrided with number of observations.
}
}
\details{
GOtree returns a GOtree object. In contains a list of significant GO terms.
plot() generated a visual plot of the GO tree.

GOtreeHits returns a vector with the genes/probes in a specific GO term.
 
GOtreeWithLeaveOut repeats function GOtree, but with different input.
GOtreeWithLeaveOut takes a table of expression data as input, performs PCA,
extracts probes / genes for the specified principal component and subsequenly
performs GOtree. This is repeated the specified number of times. It can run
with leave one out or with leave out a percentage. Only the GO terms that is
found overrepresented in all the runs "qualifies" and a new p-value is
calculated as the median of the p-value from all the runs.
An object of type GOtree is returned.
}
\value{
GOtree returns a object of type GOtree.

GOtreeWithLeaveOut returns a object of type GOtree.
}
\author{
Morten Hansen \email{mhansen@sund.ku.dk} and Jorgen Olsen
\email{jolsen@sund.ku.dk}
}
\note{
GOtreeWithLeaveOut may take some time to run - depending on the number of
samples.
}
\examples{
library(serumStimulation)
data(serumStimulation)
pcaOutput <- pca( serumStimulation )
posLoadings <- getRankedProbeIds( x=pcaOutput )
GOs <- GOtree( input=posLoadings[1:1000] )
GOs
plot(GOs, legendPosition=NULL)

\dontrun{
GOs <- GOtreeWithLeaveOut( exprsData=serumStimulation ) }
}
\keyword{ methods }

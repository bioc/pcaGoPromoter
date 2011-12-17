\name{primo}
\Rdversion{1.1}
\alias{primo}
\alias{primoHits}
\alias{primoWithLeaveOut}
\title{
primo, primoHits and primoWithLeaveOut
}
\description{
primo returns the predicted transcription factor binding sites for a given
subset of probes.

primoHits returns a list of refseqs where a given pwm scores above a threshold.

primoWithLeaveOut returns the same as \code{primo}, but run through the samples
multiple times with applied 'Leave out' cross-validation settings.
}
\usage{
primo( input, inputType = "hgu133plus2", org = "Hs",
       PvalueCutOff = 0.05, cutOff = 0.9,
       p.adjust.method = "fdr", printIgnored = FALSE , primoData = NULL )

primoHits( input, inputType = "hgu133plus2", org = "Hs",
           id, cutOff = 0.9, printIgnored = FALSE , primoData = NULL )

primoWithLeaveOut( exprsData, inputType = "hgu133plus2", org = "Hs",
                   pc = 1, decreasing = TRUE, noProbes = 1000,
                   leaveOut=1, runs = NCOL(exprsData) )
}
\arguments{
  \item{input}{
a character vector with the input. Can be Affymetrix probeset IDs, Entrez IDs
or gene symbols.
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
  \item{PvalueCutOff}{
a P-value cut off. All P-value that are below this cut off will be returned.
}
  \item{cutOff}{
the percentage of promoters that each TF should bind to in theory
}
  \item{p.adjust.method}{
the method for adjust p-values due to multiple testing. Default is "fdr"
(False Discovery Rate). See manual pages p.adjust for possibilities.
}
  \item{printIgnored}{
Should the ignored RefSeqs be printed? Default is FALSE.
}
  \item{primoData}{
primoData object from function primoData. See manual pages for primoData for
more information. Default is NULL.
}
  \item{id}{
a vector with the pwm id used to retrieve the probes with the corresponding
predicted hits. The ids should be taken from the "id" column of the primo
output.
}
  \item{exprsData}{
a table with expression data. Variables in rows and observations in columns.
Row names should be probe identifiers (Affymetrix probe set ID, Gene Symbol or
Entrez IDs), column names should be sample identifiers.
}
  \item{pc}{
a number indication which principal component to extract the loadings from.
}
  \item{decreasing}{
a logical value indication whether the loadings should be sorted in decreasing
of ascending order ( \code{decreasing == FALSE} ).
}
  \item{noProbes}{
a number indicating how many probe loadings to use
}
  \item{leaveOut}{
a number indication what percentage to leave out. If set to 1, each observation
would be left out once and runs is set equal to number of observations.
Deafault is 1.
}
  \item{runs}{
a number indicating how many times to run with leave out. If leaveOut = 1, runs
is overrided with number of observations.
}
}
\details{
primo returns the predicted transcription factor binding sites for a given set
of probes (Affymetrix probe ids, Entrez ID or gene symbols). Note: All input
values are first converted into refseqs.

primoHits returns all refseqs representing promoters where the pwm scores above
the threshold. The threshold is based on the cutOff percentage. For each
promoter the maximum value of the pwm is collected and sorted in increasing
order. The threshold is the pwm score at the cutOff percentage in the list.
Each promoter that has a scores above the threshold is considered to have
potentionel binding site of the transcription factor.

primoWithLeaveOut repeats function primo, but with different input.
primoWithLeaveOut runs primo the specified number of times. It can with leave
one out or with leave out a percentage.
Only the pwns that is represented in all the runs "qualifies" and a new p-value
is calculated as the median of the p-value from all the runs.
The list of over and under represneted pwms is returned.
}
\value{
\item{upRegulated}{
a data frame with the possible transcription factors for up regulated genes.
}
\item{downRegulated}{
a data frame with the possible transcription factors for down regulated genes.
}

	primoHits returns a character vector of probe ids (or refseqs).

	primoWithLeaveOneOut returns the same data structure as primo.
}
\references{
Stegmann, A., Hansen, M., et al. (2006). "Metabolome, transcriptome, and
bioinformatic cis-element analyses point to HNF-4 as a central regulator of gene
expression during enterocyte differentiation." Physiological Genomics 27(2):
141-155.
}
\author{
Morten Hansen \email{mhansen@sund.ku.dk} and Jorgen Olsen
\email{jolsen@sund.ku.dk}
}
\note{
The P-values in column \sQuote{pValue} of the output is not corrected for
multiple testing.

primoWithLeaveOneOut takes some time to run - especially if the is a lot of
samples.
}
\examples{
library(serumStimulation)
data(serumStimulation)
pcaOutput <- pca( serumStimulation )
negLoadings <- getRankedProbeIds( pcaOutput, pc=1, decreasing=FALSE )
TFs <- primo( negLoadings[1:1000] )
probeIds <- primoHits( negLoadings , id="M01031" )

\dontrun{
TFs <- primoWithLeaveOut( serumStimulation )
}
}
\keyword{ methods }

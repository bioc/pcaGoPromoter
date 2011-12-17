\name{primoData}
\Rdversion{1.1}
\alias{primoData}
\title{
primoData
}
\description{
primoData calculated the data objects used by primo().
}
\usage{
primoData(promoters, matrices, cleanUpPromoters = TRUE)
}
\arguments{
  \item{promoters}{
a list with promoter data. See details.
}
  \item{matrices}{
a list with matrices. See details.
}
  \item{cleanUpPromoters}{
Should some house-keeping procedures (remove dublets, map common name etc.)
be preformed on promoters? 
}
}
\details{
primoData makes a custom data object for function primo().

primoData needs 2 arguments: promoters and matrices.

promoter must consist of 2 list: a list with the RefSeq names for all the
promoters and a list with all the promoter sequences in lowercase letters.
The length of the two lists must be the same.

The function pcaGoPromoters:::primoData.getPromoter( filename ) can be used to
load promoter from a file in FASTA format.

matrices must consist of one list of lists. Each list element must have a
unique name and holds 3 data elements: baseId, name and pwm. baseId is a
character vector with a base id. name is a character vector with the common
name for the matrix. pwm is a position weighted matrix with the base A,C,G,T in
rows and the weights in columns.

Example of promoter and matrices data types is shown in the wiki at: ???

primoData returns an object of type primoData. It contains to list - threshold
and promotersRefseqs.

thresholds is a list of lists. Each list element have a unique name and holds 4
data elements: baseId, name, pwmLength and maxScores. baseId is a charactor
vector with a base id. name is a character vector with the common name for the
matrix. pwmLength is the length of the position weighted matrix. maxScores
holds the max score from all the promoters.

promotersRefseq is a list with the Refseq names for all the promoters. Some
promoters can have one that one refseq name. These names maps to the maxScores
from thresholds.
}
\value{
primoData returns an object to be used by primo() as argument primoData.
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
primoData takes very long time (many hours) to run. Install package 'multicore'
and use a computer with multiple cores.
}
\examples{
\dontrun{
myPrimoData <- primoData( myPromoters, myMatrices)
}
}
\keyword{ methods }

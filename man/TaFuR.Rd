\name{TaFuR}
\alias{TaFuR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Factor effects on microbiome community subsets contributing metagenome functions
}
\description{
TaFuR (Taxonomic and Functional Relationships) is used to characterize the microbiome community composition contributing functions that are observed or inferred across metagenomic samples.  Using CatFun the results of TaFuR are compared to results of univariate screens for differences in functions abundances across a factor of interest (e.g. Love et al. 2014) to categorize microbiomes into functional community categories.
}
\usage{
TaFuR(dat, met, no.processors = 1, method = "bray", fact = "Genus", func_col = "Gene", samp_col = "Sample", otu_col = "OTU.ID", count_col = "CountContributedByOTU")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A data frame containinng information about function by OTU by sample associations.  Column names of dat include func_col, samp_col, otu_col, and count_col.
}
  \item{met}{
A data.frame of metadata including columns samp_col and fact.
}
  \item{no.processors}{
The number of processors to use for computing individual distance matrices.
}
  \item{method}{
Any applicable distance method for vegan::vegdist.
}
  \item{fact}{
The factor of interest. Explanitory variable tested with vegan::adonis.
}
  \item{func_col}{
Column name in dat containing the gene/function names.
}
  \item{samp_col}{
Column name in dat and met containing the sample names.
}
  \item{otu_col}{
Column name in dat containing the OTU identifiers.
}
  \item{count_col}{
Column name in dat containing the counts of functions contributed to each sample by each OTU.
}
}
\details{
TaFuR parses dat into separate community matrices, one for each observed or inferred function.  Each matrix is checked to be multivariate (contributed by at least two OTUs) and that at least two levels of the factor of interest are included.   To alleviate compositionality bias the Hellinger transformation (Legendre and Gallagher 2001) is applied to individual matrices.  Distance matrices are then computed for each community subset and the effect of factor of interest on each matrix is assessed using vegan::adonis.
}
\value{
Output is a list of five with the following:
\item{1:}{per function test statistics}
\item{2:}{per function community matrices}
\item{3:}{per function adonis tables}
\item{4:}{functions private to each factor level}
\item{5:}{OTU by function matrix (this object plus data(met) can be used for DESeq2 analysis to produce data(dds))}
}

\references{
Legendre P, Gallagher ED (2001) Ecologically meaningful transformations for ordination of species data, \emph{Oecologia}. 129:271-280.

Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2, \emph{Genome Biology}, 15:550.

Phillips CD, Hanson J, Wilkinson J, Koenig L, Rees E, Webala P, Kingston T (submitted) Microbiome structural and functional incongruence explained along host dietary niche space.
}
\author{
Caleb D. Phillips <caleb.phillips@ttu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(func_dat)
data(met)
a = TaFuR(func_dat,met)
}

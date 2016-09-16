\name{CatFun}
\alias{CatFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CatFun
Frequency and community composition based catagorization of metagenome functions
}
\description{
CatFun takes the results of DESeq (Love et al. 2014) and of TaFuR {FunkyTax} to categorize functions as conserved, equivalent, divergent, or enhanced following Phillips et al. (submitted).  These categories describe how differences in metagenome functional repetoire between two or more factor levels relate to differences in underlying microbiome community structure.
}
\usage{
CatFun(tf, ds, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tf}{
Output of TaFuR, which is a list of length five.
}
  \item{ds}{
Results of DESeq2::results.
}
  \item{alpha}{
P-value threshold to be applied to Benjamini-Hochberg (Benjamini and Hochberg 1956) adjusted p-values for both DESeq2 and TAFuR analysis.  This threshold determines boundaries for function classification.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Output is a list of two with the following:
\item{1:}{ggplot object summarizing classification.}
\item{2:}{Data frame with each function rank-ordered according to abundance and with functional classifications provided.}

}
\references{
Benjamini Y, Hochberg, Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing, \emph{Journal of the Royal Statistical Society, Series B}. 57:289–300.

Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2, \emph{Genome Biology}, 15:550.

Phillips CD, Hanson J, Wilkinson J, Koenig L, Rees E, Webala P, Kingston T (submitted) Microbiome Structural and Functional Incongruence Explained along Host Dietary Niche Space.
~
}
\author{
Caleb D. Phillips
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(dds)
data(func_dat)
data(met)
a = TaFuR(func_dat, met)
b = CatFun(a, dds)

#renders plot

b[[1]]
}

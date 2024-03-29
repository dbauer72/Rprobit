% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_cl.R
\name{data_cl}
\alias{data_cl}
\title{R6 code for setting up a  'data'-object}
\value{
(valid) \code{data_cl}-object
}
\description{
\code{data_cl}-object.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{data}}{list containing the data in regressor form}

\item{\code{class_member}}{list  membership for latent classes}

\item{\code{ordered}}{Boolean indicating if choices are ordered.}

\item{\code{vars}}{string array containing the names of the regressor variables.}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{N}}{number of deciders}

\item{\code{Tp}}{choice occasions per decider}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-data_cl-new}{\code{data_cl$new()}}
\item \href{#method-data_cl-set_data}{\code{data_cl$set_data()}}
\item \href{#method-data_cl-set_class_member}{\code{data_cl$set_class_member()}}
\item \href{#method-data_cl-set_vars}{\code{data_cl$set_vars()}}
\item \href{#method-data_cl-set_ordered}{\code{data_cl$set_ordered()}}
\item \href{#method-data_cl-print}{\code{data_cl$print()}}
\item \href{#method-data_cl-plot}{\code{data_cl$plot()}}
\item \href{#method-data_cl-clone}{\code{data_cl$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-data_cl-new"></a>}}
\if{latex}{\out{\hypertarget{method-data_cl-new}{}}}
\subsection{Method \code{new()}}{
initialization function
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{data_cl$new(data = list(), ordered = FALSE, vars = "")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{list}

\item{\code{ordered}}{Boolean indicating if choices are ordered.}

\item{\code{vars}}{strings of regressor names}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-data_cl-set_data"></a>}}
\if{latex}{\out{\hypertarget{method-data_cl-set_data}{}}}
\subsection{Method \code{set_data()}}{
Set data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{data_cl$set_data(val)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{val}}{data frame}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-data_cl-set_class_member"></a>}}
\if{latex}{\out{\hypertarget{method-data_cl-set_class_member}{}}}
\subsection{Method \code{set_class_member()}}{
Set class_member
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{data_cl$set_class_member(val)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{val}}{list}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-data_cl-set_vars"></a>}}
\if{latex}{\out{\hypertarget{method-data_cl-set_vars}{}}}
\subsection{Method \code{set_vars()}}{
Set vars
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{data_cl$set_vars(val)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{val}}{vars}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-data_cl-set_ordered"></a>}}
\if{latex}{\out{\hypertarget{method-data_cl-set_ordered}{}}}
\subsection{Method \code{set_ordered()}}{
Set ordered
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{data_cl$set_ordered(val)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{val}}{ordered}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-data_cl-print"></a>}}
\if{latex}{\out{\hypertarget{method-data_cl-print}{}}}
\subsection{Method \code{print()}}{
print object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{data_cl$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-data_cl-plot"></a>}}
\if{latex}{\out{\hypertarget{method-data_cl-plot}{}}}
\subsection{Method \code{plot()}}{
plot the data frame: density plot for each chosen alternative.
plot the data set
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{data_cl$plot()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-data_cl-clone"></a>}}
\if{latex}{\out{\hypertarget{method-data_cl-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{data_cl$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

1.  Uniform formatting of `@title`.  Should be: first word only capitalized, ends with ".".
2.  `@return` lists.  These should all be of the form
```{r}
#' @return something something a list with elements: (don't forget the ":")
#' \describe{
#'   \item{\code{a}}{description of element a.}
#'   \item{\code{b}}{description of element b.}
#' }
```

3.  Latex code.  Typically this doesn't display properly in vim/emacs & rstudio.  So the trick is to define such commands twice, , e.g.:

```{r}
#' \deqn{
#' \mathbf{\Lambda}(\mathbf{Y}_t)^{\alpha - 1}
#' }{
#' \Lambda(Y_t)^(\alpha-1)
#' }
```
To see the different outputs, compile the PDF doc (this can be done from the terminal using `R CMD check LMN_XYZ.tar.gz`), and also look at output in Rstudio (which is a bit nicer than vim/emacs, but still very basic and it's what most people will be using)

4.  `@description`.  These can sometimes be made more informative, if the `@title` is not enough.

5.  Short examples for all exported functions (i.e., having the `@export` keyword).  These go in a block called `@examples`.  Also, the main doc for the package (in the file `LMN.R`) should have an example.

6.  Package vignette.  A basic walkthrough, typically using a specific example.  I'll get this started.

.onAttach <- function(lib, pkg)  {
  packageStartupMessage("biomonitoR ",
                        utils::packageDescription("biomonitoR",
                                                  fields="Version"),
                        appendLF = TRUE)
}

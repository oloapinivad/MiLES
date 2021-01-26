# cdo interface imported from ESMvaltool
cdo <-
  function(command,
           args = "",
           input = "",
           options = "",
           output = "",
           stdout = "",
           noout = F) {
    if (args != "") {
      args <- paste0(",", args)
    }
    if (stdout != "") {
      stdout <- paste0(" > '", stdout, "'")
      noout <- T
    }
    if (input[1] != "") {
      for (i in 1:length(input)) {
        input[i] <- paste0("'", input[i], "'")
      }
      input <- paste(input, collapse = " ")
    }
    output0 <- output
    if (output != "") {
      output <- paste0("'", output, "'")
    } else if (!noout) {
      output <- tempfile()
      output0 <- output
    }
    argstr <- paste0(
      options, " ", command, args, " ", input, " ", output,
      " ", stdout
    )
    print(paste("cdo", argstr))
    ret <- system2("cdo", args = argstr)
    if (ret != 0) {
      stop(paste("Failed (", ret, "): cdo", argstr))
    }
    return(output0)
  }

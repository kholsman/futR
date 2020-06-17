gdbsource2<-
function (file, interactive = FALSE) 
{
    if (!file.exists(file)) 
        stop("File '", file, "' not found")
    if (.Platform$OS.type == "windows") {
        return(.gdbsource.win(file, interactive))
        gdbscript <- tempfile()
            
       if (interactive) {
            gdbcmd <- c(paste("run --vanilla <", file), "bt")
            gdbcmd <- paste(gdbcmd, "\n", collapse = "")
            cat(gdbcmd, file = gdbscript)
            cmd <- paste("R -d gdb --debugger-args=\"-x", gdbscript, 
                "\"")
            system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = TRUE)
            return(NULL)
        }
        else {
            cat("run\nbt\nquit\n", file = gdbscript)
            cmd <- paste("R --vanilla < ", file, " -d gdb --debugger-args=\"-x", 
                gdbscript, "\"")
            txt <- system(cmd, intern = TRUE, ignore.stdout = FALSE, 
                ignore.stderr = TRUE)
            attr(txt, "file") <- file
            class(txt) <- "backtrace"
            return(txt)
        }

    }else{

        gdbscript <- tempfile()
        
        if (interactive) {
            gdbcmd <- c(paste("run --vanilla <", file), "bt")
            gdbcmd <- paste(gdbcmd, "\n", collapse = "")
            cat(gdbcmd, file = gdbscript)
            cmd <- paste("R -d lldb --debugger-args=\"-x", gdbscript, 
                "\"")
            system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = TRUE)
            return(NULL)
        }
        else {
            cat("run\nbt\nquit\n", file = gdbscript)
            cmd <- paste("R --vanilla < ", file, " -d lldb --debugger-args=\"-x", 
                gdbscript, "\"")
              
            txt <- system(cmd, intern = TRUE, ignore.stdout = FALSE, 
                ignore.stderr = TRUE)
             attr(txt, "file") <- file
             #class(txt) <- "backtrace"
            return(txt)
        }
    }
}
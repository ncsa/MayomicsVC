###########################################################################################


###########################################################################################

# The Task block is where the variables and the functions are defined for performing a certain task

task HelloTask {
        
   String Name		   # Name of the Sample
   String BashScriptPath

   command <<<
 
      /bin/bash ${BashScriptPath} ${Name}
 
    >>>

   output {

   File out = stdout()

   }
   

}  # End of task block


workflow OutputFlow {

   call HelloTask

}

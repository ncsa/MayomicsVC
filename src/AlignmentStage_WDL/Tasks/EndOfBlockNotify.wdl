#######################################################################################
##      This wdl script sends an email to analyst at the end of each step      ##
#######################################################################################

task EndOfStageEmailTask {

   Array[Int] Email
   String Failure_Logs
   String dollar = "$"


   command <<<
     
      flag=${Email[0]}

      if [ ${dollar}{flag} -eq 0 ]
      then
         echo "The Alignment Block has completed. The attached Failure_logs.txt has information on individual samples" | mailx -s "End of Alignment Block" rvenka21@illinois.edu
      ##elif[${sep=',' Email} -eq 1 ]
      ##then
        ## echo "The ReCal Block has completed. The attached Failure_logs.txt has information on individual samples" | mailx -a ${Failure_Logs} -s "End of ReCal Block" rvenka21@illinois.edu 

      ##elif[${sep=',' Email} -eq 2 ]
      ##then
        ## echo "The Variant Calling Block has completed. The attached Failure_logs.txt has information on individual samples" | mailx -a ${Failure_Logs} -s "End of Variant Calling block" rvenka21@illinois.edu
      fi
   >>>

   runtime {
      continueOnReturnCode: false
   }

}

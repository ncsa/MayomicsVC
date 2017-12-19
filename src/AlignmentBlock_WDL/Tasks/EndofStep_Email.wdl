#######################################################################################
##      This wdl script sends an email to analyst at the end of each step      ##
#######################################################################################

task EndofStep_Email {

   Array[Int] Email
   String Failure_Logs

   command {

      if[ ${Email} -eq 0 ]
      then
         echo "The BWA Step has completed. The attached Failure_logs.txt has information on individual samples" | mailx -a ${Failure_Logs} -s "End of BWA Step" rvenka21@illinois.edu 

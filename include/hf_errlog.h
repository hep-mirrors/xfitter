#pragma once
#include<string>
/// Log errors, warnings and info messages
/*
The severity level is coded in the first two characters of the string message:
"I:" - information   (continue execution always)
"W:" - warning       (continue execution by default)
"S:" - serious error (terminate program  by default)
"F:" - fatal error   (terminate program  always)

id is an error number
To add a new error number so that it would not clash with already 
existing ones it is proposed to put the date of coding into it:        
               error number = YYMMDDnn                            
         (e.g. 120117nn for the date = 17.01.2012)                 
This gives for code developers 100 numbers per day which is enough.
You will never clash with yourself and quite unlikely with others.
For imported code, reserve 100 IDs for every imported module 
*/
void hf_errlog(int id,const std::string&message);

#pragma once
#include<string>
/// Report errors and messages with different level of severity (I:, W:, S:). S: stops the execution
void hf_errlog(int id,const std::string&message);

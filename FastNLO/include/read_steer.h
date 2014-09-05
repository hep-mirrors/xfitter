// Author: Daniel Britzger
// DESY, 17/07/2012
#ifndef __READ_STEER_INC__
#define __READ_STEER_INC__

/**
// **********************************************************************************
//
//     read_steer.h
//     Tiny reading tool to read steering values from one or more steering files.
//
//     This class reads in values, which are stored in a file. New variables
//     can be included without changes of the steering class.
//
//     Features
//     ------------------------------
//       o  Following types are supported:
//            - Single values
//              bool, int, double, string (with empty spaces), char
//            - Arrays
//              int, double, string (with empty spaces)
//            - Tables/matrices
//              int, double, string
//       o  Multiple files can be read in and handled individually or together.
//       o  Namespaces can be defined (e.g. same variable name in different namesapces).
//       o  Variables within a steer-file can be defined similar to shell skripts.
//       o  Command line arguments can be parsed and can superseed values from files.
//       o  Easy access via pre-processor commands
//       o  Other files can be included into steering file
//
//
//     Initalize the steering
//     ------------------------------
//     Set the filename and initilize the read_steer class by using:
//        read_steer::readfile(string filename);
//     or
//        READ("steerfile.str");
//
//
//     Single values
//     ------------------------------
//     To access the values, simply use the adequate getter functions
//     for the desired variable type and the label of this variable.
//     To speed up the code and to avoid repeated string comparisions,
//     use static variables, e.g if you want to access the value in
//     your steering file with the label 'pi' or 'name', use:
//        static double pi   = read_steer::getdouble("pi");
//     or
//        static string name = read_steer::getstring("name");
//        static int    age  = read_steer::getint("age");
//        static bool   sex  = read_steer::getbool("female");
//     A more readable method is to use the pre-processor macro
//        static string name = STRING(name);
//        static int     age = INT(age);
//        static bool    sex = BOOL(female);
//
//     Labels are case sensitive!
//
//
//
//     Syntax of steering file
//     ------------------------------
//     The steering file can consist of an arbitrary number of lines, where
//     the syntax should follow:
//         <label>              <value>         [!comment]
//     where 'label' and 'value' are necessary tags and comments are
//     beginning with the '!' character and are ignored by the read_steer class.
//     As seperator between the <label> and the <value> empty spaces or tabstops
//     are recognized. Complete lines can beginn with '!' to mark comments.
//     If string values should contain empty spaces, enclose them in double quotes
//     like:
//         Name                 "Peter Higgs"
//         Age                  137
//     Boolean values can be assigned by 0, 1, true or false, e.g.
//         WithHiggs            true
//
//
//     Arrays
//     ------------------------------
//     To read in an array of values, assign a label and enclose the following
//     values in curly brackets { }, with a leading empty space [" {" and "}"].
//     Within curly brackets, each separated (by an empty space or tabstop)
//     value is stored in the array as an element, or each line when double quotes are
//     used (only one occurence of double quotes per line is recognized).
//     The steering file should look like:
//
//        !Numbers 1-9 are read into an array (9 elements)
//        Array1 {                      !array starts here
//          1 2 3 4 5                   ! integers from 1 to 5
//          6 7 8 9
//        }
//        !Eleven Names of famous musicians (11 elements)
//        FamousMusicians {
//          John Paul Ringo George      !The Beatles
//          Beethoven Bach Mozart       !Famous componists
//          Mick Keith Ron Charlie      !The Rolling Stones
//        }
//        !Full sentences or documentations (2 string-elements)
//        Array3 {
//          "Hello World!"              !Sentences are great
//          "Was the first scream."
//        }
//
//     To access the arrays use e.g.:
//         static vector<string> musicians = read_steer::getstringarray("FamousMusicians");
//         static vector<double> nums      = read_steer::getdoublearray("Array1");
//         static vector<int>    ints      = read_steer::getintarray("Array1");
//     or equivalently
//         static vector<string> musicians = STRING_ARR(FamousMusicians);
//         static vector<double> nums      = DOUBLE_ARR(Array1);
//         static vector<int>    ints      = INT_ARR(Array1);
//
//
//     Tables and matrices
//     ------------------------------
//     Tables and matrices are tagged by ' {{' and '}}' in the steering file.
//     No double quotes (e.g. "text") are allowed as table values. The first row
//     of a table is always expected to be the row-header. The row headers are
//     separated by empty spaces or tabstops. If matrices are necessary
//     keep the first line empty or add a comment there. The steering file should look like:
//
//        Crossections {{
//              Q2min   Q2max   cs[pb]  stat[%]                 !header tags should not contain emtpy spaces
//              100     200     22.12   1.2
//              200     300     12.72   2.7
//              300     500     23.22   5.3
//        }}
//        Participants {{
//              Name            Surname         Country         ! first line is always the header
//              Obama           Barack          U.S.A.
//              Merkel          Angela          Germany
//              Benedikt        XVI             Vatican
//        }}
//        Matrix {{
//              !the first line is ignored. Keep it empty.
//              11 12
//              21 22
//         }}
//
//      To access the table use e.g.:
//          static vector<vector<string> > guys = read_steer::getstringtable("Participants");
//          static vector<vector<double> > cs   = read_steer::getdoubletable("Crossections");
//          static vector<vector<int> >    mat  = read_steer::getinttable("Matrix");
//      To access the table header use:
//         static vector<string>           head = read_steer::gettableheader("Crossections);
//      To access a single column of a table use:
//          static vector<double> xs            = read_steer::getdoublecolumn("Crossections","cs[pb]");
//          static vector<string> nick          = read_steer::getstringcolumn("Participants","Surname");
//     or equivalently
//          static vector<double> xs            = DOUBLE_COL("Crossections","cs[pb]");
//          static vector<string> nick          = STRING_COL("Participants","Surname");
//
//
//
//     Multiple steering files.
//     ------------------------------
//     In case multiple steering files are necessary, each steering file must
//     be assigned a unique 'steerID' if variable names (labels)  are identically.
//          read_steer::readfile("file1.steer","file1")
//          read_steer::readfile("anotherfile.steer","constants")
//
//     To access values, pass the steerID to the getter methods, e.g.:
//          static double pi   = read_steer::getdouble("pi","constants");
//          static string name = read_steer::getstring("name","file1");
//          static vector<vector<string> > ConfIchepNames = read_steer::getstringcolumn("Participants","Surname","file1")
//     You can access the values at any place within your code.
//
//     If different labels should be read in from multiple files, just call
//          read_steer::readfile("file1.steer");
//          read_steer::readfile("file2.steer");
//     and access the variables without the usage of the steerID.
//
//
//     Namespaces
//     ------------------------------
//     Instead of using multiple files for reading identical labels for
//     various occasions, one can use namespaces instead. Namespaces are
//     handled identically to multiple files, but can be defined within
//     one single steering file. Each namespace is assigned a steerID.
//     Namespaces are defined by a label, which is used as the steerID
//     and start with the '{{{' tag and end with the '}}}' tag.
//     A steerfile could look like:
//       HostInstitute          CERN            ! standard variabel
//       ATLAS {{{                              ! namespace ATLAS starts here
//          length      45                      ! define variables as usual
//          height      22
//          weight      7000
//          Crossection {{                      ! also tables are possible
//              bin     cs[pb]  stat[%]
//              1       32.2    1.2
//              2       12.2    3.2
//          }}
//       }}}                                    ! namespace ATLAS ends here
//       CMS {{{
//          length      21
//          height      16
//          weight      12500
//          Crossection {{
//              bin     cs[pb]  stat[%]
//              1       33.1    0.8
//              2       13.6    3.4
//          }}
//       }}}
//
//     To access the values, use the steerID which is the label of the namespace
//          static double ATLASheight     = read_steer::getdouble("height","ATLAS");
//          static double CMSheight       = read_steer::getdouble("height","CMS");
//          static vector<double> CMSxs   = read_steer::getdoublecolumn("Crossection","cs[pb]","CMS");
//          static vector<double> ATLASxs = read_steer::getdoublecolumn("Crossection","cs[pb]","ATLAS");
//
//     Warning: Namespace steerID and file steerID might conflict if identically!
//     It is NOT possible to read in multiple files, wherein identical namespaces are define!
//     There is no possiblity to access variables using the pre-processor commands (e.g. DOUBLE(val))
//     for namespaces, other than the standard namespace.
//
//
//     Script-like Variables
//     ------------------------------
//     It is often neessary to read in identical substrings, e.g. if
//     many different files are located in the same folder. To simplify the
//     structure of the steering file, 'script-like' variables can be used.
//
//     It is possible to access previously defined variables foo by ${foo}.
//     An example steering file can look like:
//           !Home directories of famous physicists
//           HomeDir                    /afs/cern.ch/user
//           UserEinstein               ${HomeDir}/e/einstein
//           UserNewton                 "${HomeDir}/i/isaac"
//
//     Local variables are only valid within the defined namespace.
//
//
//     Parse command line
//     ------------------------------
//     It is possible to read in steering values and specify steering files
//     via the command line of the program.
//     To read in values over the command line, one has to call
//            read_steer::parsecommandline(argc,argv);
//        or
//            PARSE(arc,argv);
//     Where argc and argv are the command line parameters as specified in main(argc,argv).
//     Specify a value when executing the program (e.g. Run) over the command line like:
//            >$ Run label1=value1 WelcomeMessage="Hello World" Name=Einstein:Names steerfile=file.str steerfile=file2.str:Names
//
//     This example will initialize the labels in the standard namespace
//            label1               value1
//            WelcomeMessage       "Hello World"
//     Further it will read in the steering file file.str into the standard namespace
//     and will read in the file file2.str into the namespace "Names".
//     The value
//            Name                 Einstein
//     will be available in the namespace 'Names'.
//
//
//     Include external files into steer-file
//     ---------------------------------------
//     It is possible to include other files into a steerfile. This might be useful
//     for defining tables or arrays.
//     To include other files, use the '#include:<filename>[:start[:stop]]' tag like:
//            #include:steerfile2.str
//     The file steerfile2.str in this example is handled like all its content would be at
//     this position of the base-steerfile. This is useful if e.g. cross section
//     tables should be read in:
//            CrossSections {{
//                #include:HiggsCrossSection.txt
//            }}
//     If only certain lines of an external file should be read in, these could be
//     specified, with a separated ':'-symbol.
//            CrossSections {{
//                #include:HiggsCrossSection.txt:2:12
//            }}
//      Here, only lines 2-12 are read in. It is also possible, to specify only the first line.
//
//      Mention: Although the content of the included file is exactly treated like its content
//      would stand in the base steefile, it is not possible to define namespaces within included files.
//      However, brackets, local variables, table definitions, etc. are all treated the same way.
//
//
//
//     Printing
//     ------------------------------
//     Print all steering information in SingleFileMode by calling
//          read_steer::print();
//     If multiple files are used, print all information using:
//          read_steer::printall();
//        or
//          PRINTALL();
//     or just the information of one steerID:
//          read_steer::print("constants");
//
//
//     Warning
//     ------------------------------
//     This tool is basically based on string comparisions to identify
//     the steering values. Do not call it too often within your code in order to
//     avoid speed problems. The best way to use the steering values is to assign the
//     values to static const variables.
//
//
//     Not existing labels or steerID
//     ------------------------------
//     If you access an element which was not read in from the steering file,
//     this label is automatically added to the list of elements with value zero.
//     If a steerID is accessed, which is not identified, a new steerID is added
//     to the list without any values.
//
//
//     D. Britzger
//     daniel.britzger@desy.de
//
//
// ********************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
//update to unordered_map in C++11 in gcc4.7
#include <vector>

// use the pre-processor for accessing values
// in 'single-file' mode
#define READ(X) read_steer::readfile(X)

#define BOOL(X) read_steer::getbool(#X)
#define INT(X) read_steer::getint(#X)
#define DOUBLE(X) read_steer::getdouble(#X)
#define CHAR(X) read_steer::getstring(#X).c_str()
#define STRING(X) read_steer::getstring(#X)

#define BOOL_ARR(X) read_steer::getboolarray(#X)
#define INT_ARR(X) read_steer::getintarray(#X)
#define DOUBLE_ARR(X) read_steer::getdoublearray(#X)
#define STRING_ARR(X) read_steer::getstringarray(#X)

#define BOOL_COL(X,Y) read_steer::getboolcolumn(#X,#Y)
#define INT_COL(X,Y) read_steer::getintcolumn(#X,#Y)
#define DOUBLE_COL(X,Y) read_steer::getdoublecolumn(#X,#Y)
#define STRING_COL(X,Y) read_steer::getstringcolumn(#X,#Y)

#define INT_TAB(X) read_steer::getinttable(#X)
#define DOUBLE_TAB(X) read_steer::getdoubletable(#X)
#define STRING_TAB(X) read_steer::getstringtable(#X)

// use the pre-processor for accessing values
// in 'multi-file' mode
// NS specifies the steerID-namespace
#define READ_NS(X,NS) read_steer::readfile(X,NS)

#define BOOL_NS(X,NS) read_steer::getbool(#X,NS)
#define INT_NS(X,NS) read_steer::getint(#X,NS)
#define DOUBLE_NS(X,NS) read_steer::getdouble(#X,NS)
#define CHAR_NS(X,NS) read_steer::getstring(#X,NS).c_str()
#define STRING_NS(X,NS) read_steer::getstring(#X,NS)

#define BOOL_ARR_NS(X,NS) read_steer::getboolarray(#X,NS)
#define INT_ARR_NS(X,NS) read_steer::getintarray(#X,NS)
#define DOUBLE_ARR_NS(X,NS) read_steer::getdoublearray(#X,NS)
#define STRING_ARR_NS(X,NS) read_steer::getstringarray(#X,NS)

#define BOOL_COL_NS(X,Y,NS) read_steer::getboolcolumn(#X,#Y,NS)
#define INT_COL_NS(X,Y,NS) read_steer::getintcolumn(#X,#Y,NS)
#define DOUBLE_COL_NS(X,Y,NS) read_steer::getdoublecolumn(#X,#Y,NS)
#define STRING_COL_NS(X,Y,NS) read_steer::getstringcolumn(#X,#Y,NS)

#define INT_TAB_NS(X,NS) read_steer::getinttable(#X,NS)
#define DOUBLE_TAB_NS(X,NS) read_steer::getdoubletable(#X,NS)
#define STRING_TAB_NS(X,NS) read_steer::getstringtable(#X,NS)

// shortcuts for parsing command line and printing
#define PARSE(X,Y) read_steer::parsecommandline(X,Y)
#define PRINTALL() read_steer::printall()



using namespace std;

class read_steer {

public:
   static const string stdID;

public:
   ~read_steer() {;};
   static read_steer* Steering(string steerID=read_steer::stdID);                       // get an object!
   static void destroy();                                               // destroy all instances

public:
   // getters for single instance
   // values
   bool getb(string label);
   int geti(string label);
   double getd(string label);
   string gets(string label);
   // arrays
   vector<bool> getbf(string label);
   vector<int> getif(string label);
   vector<double> getdf(string label);
   vector<string> getsf(string label);
   // tables/matrices
   vector<string> getsthead(string label);
   vector<vector<int> > getit(string label);
   vector<vector<double> > getdt(string label);
   vector<vector<string> > getst(string label);
   vector<bool> getbtcol(string label,string col);
   vector<int> getitcol(string label,string col);
   vector<double> getdtcol(string label,string col);
   vector<string> getstcol(string label,string col);

   // setter
   void addlabel(const string lab, const string val);
   // controls
   void inits(string filename);
   int initnmspc(ifstream& strm, string filename);
   void prt();
   static void initnamespace(ifstream& strm,string filename, string steerID=read_steer::stdID) {        // set the steer-filename
      read_steer::Steering(steerID)->initnmspc(strm,filename); }


private:
   read_steer();
   read_steer(const read_steer& ) {;};
   static map<string,read_steer*>* instances;

   int read_stdin(string filename);
   int readstrm(ifstream& strm,unsigned int lstart=0,unsigned int lend=0, bool incfile=false);
   bool ParseString(string value);
   bool ParseFindString(const string str, const string tag) const;
   string ParseEnclosedString(const string) const;
   bool EnclosedStringToOneEntity(string&) const;
   int ReplaceVariables(string& value);
   bool CheckNumber(const string str) const;
   bool CheckInt(const string str) const;
   bool StringToBool(const string str, const string label="") const;
   static int cmdlinetag(const char* arg, string& label, string& value);
   static int separatetag(string& vallhs, string& valrhs, const string sep);

   map<string,string> fstrings;
   map<string,vector<string> > ffields;
   map<string,vector<vector<string> > > ftables;
   map<string,vector<string> > ftableheaders;
   bool fParseFieldMode;
   int  fParseTableMode;
   string ffieldlabel;
   vector<string> ffieldvalues;
   vector<vector<string> > ftablevalues;
   string ffilename;
   string fcurrentfilename;
   ifstream ffile;

   const string str_sep;
   const string str_cmt;
   const string str_arrbeg;
   const string str_arrend;
   const string str_tabbeg;
   const string str_tabend;
   const string str_nmspcbeg;
   const string str_nmspcend;
   const string str_inc;
   int fParseIncMode;
   const string oW;
   const string oI;
   const string oE;

public:
   // static member function
   static void readfile(string filename,string steerID=read_steer::stdID) {     // set the steer-filename
      read_steer::Steering(steerID)->inits(filename); }
   // getters
   // values
   static bool getbool(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getb(label); }
   static int getint(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->geti(label); }
   static double getdouble(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getd(label); }
   static string getstring(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->gets(label); }
   // arrays
   static vector<bool> getboolarray(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getbf(label); }
   static vector<int> getintarray(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getif(label); }
   static vector<double> getdoublearray(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getdf(label); }
   static vector<string> getstringarray(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getsf(label); }
   // tables header
   static vector<string> gettableheader(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getsthead(label); }
   // tables/matrices
   static vector<vector<int> > getinttable(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getit(label); }
   static vector<vector<double> > getdoubletable(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getdt(label); }
   static vector<vector<string> > getstringtable(string label,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getst(label); }
   // table columns
   static vector<bool> getboolcolumn(string label,string column,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getbtcol(label,column); }
   static vector<int> getintcolumn(string label,string column,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getitcol(label,column); }
   static vector<double> getdoublecolumn(string label,string column,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getdtcol(label,column); }
   static vector<string> getstringcolumn(string label,string column,string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getstcol(label,column); }

   static void printall();                                              // print values of all files
   static void print(string steerID=read_steer::stdID);                         // print values
   static bool parsecommandline(int argc,char** argv);


};

#endif

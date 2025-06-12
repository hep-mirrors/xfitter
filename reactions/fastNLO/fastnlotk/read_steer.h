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
//     Labels are case sensitive.
//
//     Check existence of label using:
//         bool IsPresent = EXIST(name);
//         bool IsPresent = read_steer::exist("name");
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
//     be assigned a unique 'steerID' if variable names (labels) are identical.
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
//     Check existence of labels in namespaces using:
//         bool IsPresent = EXIST_NS(name,"ATLAS");
//         bool IsPresent = read_steer::getexist("name","ATLAS");
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
//            read_steer::parsecommandline(vector<string> v);
//        or
//            PARSE(arc,argv);
//            PARSEV(v);
//     Where argc and argv are the command line parameters as specified in main(argc,argv).
//     Specify a value when executing the program (e.g. Run) over the command line like:
//            >$ Run label1=value1 WelcomeMessage="Hello World" Names::Name=Einstein steerfile=file.str steerfile=file2.str->Names
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
//      Adding values during runtime with setters
//      -------------------------------------------
//      It is possible to add values with a key during runtime, i.e. using functions
//      instead of a steering file. Using this option it is also possible
//      to modify and replace values.
//      For adding a single value, array or table, the shorthand notations are respectively:
//         ADD(K,Y);
//         ADDARRAY(K,Y);
//         ADDTABLE(K,H,Y);
//      For instance to add a value "76" for the key "age"
//         ADD("age",76);
//      The full notation woule be:
//         read_steer::Steering()->AddLabel("age",76);
//
//      To add a key to a specific namespace use ADD_NS, ADDARRAY_NS or ADDTABLE_NS. For instance:
//         ADD_NS("age",76,"MyMom")
//      which adds the value 76 for the label 'age' to the namespace 'MyMom'.
//
//      Further methods for appending single elements to arrays or tables also exist.
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
#include <set>
//update to unordered_map in C++11 in gcc4.7
#include <vector>
#include <sstream>
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
#define TABLEHEADER(X) read_steer::gettableheader(#X)

#define INT_TAB(X) read_steer::getinttable(#X)
#define DOUBLE_TAB(X) read_steer::getdoubletable(#X)
#define STRING_TAB(X) read_steer::getstringtable(#X)

#define ADD(X,Y) read_steer::addvalue(X,Y)
#define ADDARRAY(X,Y) read_steer::addarray(X,Y)
#define ADDTABLE(X,Y,Z) read_steer::addtable(X,Y,Z)

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
#define TABLEHEADER_NS(X,NS) read_steer::gettableheader(#X,NS)

#define INT_TAB_NS(X,NS) read_steer::getinttable(#X,NS)
#define DOUBLE_TAB_NS(X,NS) read_steer::getdoubletable(#X,NS)
#define STRING_TAB_NS(X,NS) read_steer::getstringtable(#X,NS)

#define ADD_NS(X,Y,NS) read_steer::addvalue(X,Y,NS)
#define ADDARRAY_NS(X,Y,NS) read_steer::addarray(X,Y,NS)
#define ADDTABLE_NS(X,Y,Z,NS) read_steer::addtable(X,Y,Z,NS)

// shortcuts for parsing command line and printing
#define PARSE(X,Y) read_steer::parsecommandline(X,Y)
#define PARSEV(X) read_steer::parsecommandline(X)
#define PRINTALL() read_steer::printall()

// check if key exists
#define EXIST(X) read_steer::getexist(#X)
#define EXIST_NS(X,NS) read_steer::getexist(#X,NS)

// check if array or array element exists
#define EXISTARRAY_NS(X,NS) read_steer::getarrayexist(#X, NS)
#define CONTAINKEYARRAY_NS(X,A,NS) read_steer::getarraycontainkey(#X, #A, NS)
#define PUSHBACKARRAY_NS(X,A,NS) read_steer::arraypushback_steer(X, #A, NS)

class read_steer {

private:
   static std::string stdID;

public:
   ~read_steer() {;};
   static read_steer* Steering(std::string steerID=read_steer::stdID);                       // get an object!
   static void destroy();                                               // destroy all instances

public:
   // getters for single instance
   // values
   bool getb(const std::string& label);
   int geti(const std::string& label);
   double getd(const std::string& label);
   std::string gets(const std::string& label);
   // arrays
   std::vector<bool> getbf(const std::string& label);
   std::vector<int> getif(const std::string& label);
   std::vector<double> getdf(const std::string& label);
   std::vector<std::string> getsf(const std::string& label);
   // tables/matrices
   std::vector<std::string> getsthead(const std::string& label);
   std::vector<std::vector<int> > getit(const std::string& label);
   std::vector<std::vector<double> > getdt(const std::string& label);
   std::vector<std::vector<std::string> > getst(const std::string& label);
   std::vector<bool> getbtcol(const std::string& label,const std::string& col);
   std::vector<int> getitcol(const std::string& label,const std::string& col);
   std::vector<double> getdtcol(const std::string& label,const std::string& col);
   std::vector<std::string> getstcol(const std::string& label,const std::string& col);

   // check if key exists
   bool exist(const std::string& label) {
      bool ret = ( fstrings.count(label) > 0 || ffields.count(label) > 0 || ftables.count(label) > 0 );
      return ret;}
   bool arrayexist(const std::string& label) {
      return ffields.count(label) > 0;}
   bool arraycontainkey(const std::string& key, const std::string& label) {
      if (ffields.count(label) == 0) return false;
      for (size_t i = 0; i < ffields[label].size(); ++i) {
        if (ffields[label][i].find(key+"=") == std::string::npos) {
          return true;
        }
      }
      return false;}
   void arraypushback(const std::string& value, const std::string& label) {
    ffields[label].push_back(value);}
   // setter
   void AddLabel(const std::string& key, const std::string& value);
   template <typename T> void AddLabel ( const std::string& key, T value);
   void AddArray(const std::string& key, const std::vector<std::string>& values);
   template <typename T> void AddArray ( const std::string& key, const std::vector<T>& values);
   void AddTable(const std::string& key, const std::vector<std::string>& header, const std::vector<std::vector<std::string> >& values);
   template <typename T> void AddTable ( const std::string& key, const std::vector<std::string>& header, const std::vector<std::vector<T> >& values);
   void AppendToArray(const std::string& key, const std::string& entry);
   template <typename T> void AppendToArray ( const std::string& key, const T& entry);
   void AppendToTable(const std::string& key, const std::vector<std::string>& entry);
   template <typename T> void AppendToTable ( const std::string& key, const std::vector<T>& entry);
   // controls
   // report return code instead of void function
   int inits(std::string filename);
   int initnmspc(std::ifstream& strm, std::string filename);
   void prt();
   static void initnamespace(std::ifstream& strm,std::string filename, std::string steerID=read_steer::stdID) {        // set the steer-filename
      read_steer::Steering(steerID)->initnmspc(strm,filename); }

   // get labels
   std::set<std::string> GetAvailableLabels() const;
   std::set<std::string> GetAvailableArrrays() const;
   std::set<std::string> GetAvailableTables() const;

   static bool CheckNumber(const std::string& str);
   static bool CheckInt(const std::string& str);
   static int separatetag(std::string& vallhs, std::string& valrhs, const std::string& sep);

   static const std::string& GetDefaultNamespace(){ return stdID;}
   static void SetDefaultNamespace(const std::string& nspc){ stdID = nspc;}

private:
   read_steer();
   read_steer(const read_steer& ) {;};
   static std::map<std::string,read_steer*>* instances;

   int read_stdin(const std::string& filename);
   int readstrm(std::ifstream& strm,unsigned int lstart=0,unsigned int lend=0, bool incfile=false);
   bool ParseString(std::string value);
   bool ParseFindString(const std::string& str, const std::string& tag) const;
   std::string ParseEnclosedString(const std::string&) const;
   bool EnclosedStringToOneEntity(std::string&) const;
   int ReplaceVariables(std::string& value);
   bool StringToBool(const std::string& str, const std::string& label="") const;
   static int cmdlinetag(const char* arg, std::string& label, std::string& value);

   std::map<std::string,std::string> fstrings;
   std::map<std::string,std::vector<std::string> > ffields;
   std::map<std::string,std::vector<std::vector<std::string> > > ftables;
   std::map<std::string,std::vector<std::string> > ftableheaders;
   bool fParseFieldMode;
   int  fParseTableMode;
   std::string ffieldlabel;
   std::vector<std::string> ffieldvalues;
   std::vector<std::vector<std::string> > ftablevalues;
   std::string ffilename;
   std::string fcurrentfilename;
   std::ifstream ffile;

   const std::string str_sep;
   const std::string str_cmt;
   const std::string str_arrbeg;
   const std::string str_arrend;
   const std::string str_tabbeg;
   const std::string str_tabend;
   const std::string str_nmspcbeg;
   const std::string str_nmspcend;
   const std::string str_inc;
   int fParseIncMode;
   const std::string oW;
   const std::string oI;
   const std::string oE;

public:
   // static member function; report return code instead of void function
   static int readfile(std::string filename,std::string steerID=read_steer::stdID) {     // set the steer-filename
      return read_steer::Steering(steerID)->inits(filename); }
   // setters
   // verbosity
   static int fVerbosity;
   static void setVerbosity(int iVerbosity) {
      read_steer::fVerbosity = iVerbosity; }
   // getters
   // verbosity
   static int getVerbosity() {
      return read_steer::fVerbosity; }
   // values
   static bool getbool(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getb(label); }
   static int getint(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->geti(label); }
   static double getdouble(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getd(label); }
   static std::string getstring(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->gets(label); }
   // arrays
   static std::vector<bool> getboolarray(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getbf(label); }
   static std::vector<int> getintarray(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getif(label); }
   static std::vector<double> getdoublearray(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getdf(label); }
   static std::vector<std::string> getstringarray(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getsf(label); }
   // tables header
   static std::vector<std::string> gettableheader(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getsthead(label); }
   // tables/matrices
   static std::vector<std::vector<int> > getinttable(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getit(label); }
   static std::vector<std::vector<double> > getdoubletable(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getdt(label); }
   static std::vector<std::vector<std::string> > getstringtable(std::string label,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getst(label); }
   // table columns
   static std::vector<bool> getboolcolumn(std::string label,std::string column,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getbtcol(label,column); }
   static std::vector<int> getintcolumn(std::string label,std::string column,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getitcol(label,column); }
   static std::vector<double> getdoublecolumn(std::string label,std::string column,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getdtcol(label,column); }
   static std::vector<std::string> getstringcolumn(std::string label,std::string column,std::string steerID=read_steer::stdID) {
      return read_steer::Steering(steerID)->getstcol(label,column); }
   // check existence of value
   static bool getexist(const std::string& label, std::string steerID=read_steer::stdID ){
      //      if ( read_steer::instances==NULL ) read_steer::instances = new std::map<string,read_steer*>();
      if ( read_steer::instances==NULL ) return false;
      if ( read_steer::instances->count(steerID) == 0 ) return false;
      return read_steer::Steering(steerID)->exist(label);}
   static bool getarrayexist(const std::string& label, std::string steerID=read_steer::stdID ){
      //      if ( read_steer::instances==NULL ) read_steer::instances = new std::map<string,read_steer*>();
      if ( read_steer::instances==NULL ) return false;
      if ( read_steer::instances->count(steerID) == 0 ) return false;
      return read_steer::Steering(steerID)->arrayexist(label);}
   static bool getarraycontainkey(const std::string& key, const std::string& label, std::string steerID=read_steer::stdID ){
      //      if ( read_steer::instances==NULL ) read_steer::instances = new std::map<string,read_steer*>();
      if ( read_steer::instances==NULL ) return false;
      if ( read_steer::instances->count(steerID) == 0 ) return false;
      return read_steer::Steering(steerID)->arraycontainkey(key, label);}
   static void arraypushback_steer(const std::string& value, const std::string& label, std::string steerID=read_steer::stdID ){
      //      if ( read_steer::instances==NULL ) read_steer::instances = new std::map<string,read_steer*>();
      if ( read_steer::instances==NULL ) return;
      if ( read_steer::instances->count(steerID) == 0 ) return;
      read_steer::Steering(steerID)->arraypushback(value, label);}
   // add values
   template <typename T>
   static void addvalue(const std::string& key,const T& val,std::string steerID=read_steer::stdID) {
      read_steer::Steering(steerID)->AddLabel(key,val); }
   template <typename T>
   static void addarray(const std::string& key,const std::vector<T>& val,std::string steerID=read_steer::stdID) {
      read_steer::Steering(steerID)->AddArray(key,val); }
   template <typename T>
   static void addtable(const std::string& key, const std::vector<std::string>& header, const std::vector<std::vector<T> >& values,std::string steerID=read_steer::stdID) {
      read_steer::Steering(steerID)->AddTable(key,header,values); }
   template <typename T>
   static void appendtoarray(const std::string& key, const T& entry, std::string steerID=read_steer::stdID) {
      read_steer::Steering(steerID)->AppendToArray(key,entry); }
   template <typename T>
   static void appendtotable(const std::string& key, const std::vector<T>& entry, std::string steerID=read_steer::stdID) {
      read_steer::Steering(steerID)->AppendToTable(key,entry); }

   static void printall();                                              // print values of all files
   static void print(std::string steerID=read_steer::stdID);                         // print values
   static bool parsecommandline(int argc,char** argv);
   static bool parsecommandline(std::vector<std::string> argv);

};


template <typename T>
void read_steer::AddLabel ( const std::string& key, T val) {
   std::stringstream ss;
   ss << val;
   AddLabel(key,ss.str());
}

template <typename T>
void read_steer::AddArray ( const std::string& key, const std::vector<T>& val) {
   std::vector<std::string> str(val.size());
   for ( unsigned int i = 0 ; i<val.size() ; i++ ) {
      std::stringstream ss;
      ss << val[i];
      str[i] = ss.str();
   }
   AddArray(key,str);
}

template <typename T>
void read_steer::AddTable ( const std::string& key, const std::vector<std::string>& header, const std::vector<std::vector<T> >& values) {
   std::vector<std::vector<std::string> > str(values.size());
   for ( unsigned int i = 0 ; i<values.size() ; i++ ) {
      str[i].resize(values[i].size());
      for ( unsigned int j = 0 ; j<values[i].size() ; j++ ) {
         std::stringstream ss;
         ss << values[i][j];
         str[i][j]=ss.str();
      }
   }
   AddTable(key,header,str);
}

template <typename T>
void read_steer::AppendToArray ( const std::string& key, const T& entry ) {
   std::stringstream ss;
   ss << entry;
   std::string str=ss.str();
   AppendToArray(key,entry);
}

template <typename T>
void read_steer::AppendToTable ( const std::string& key, const std::vector<T>& entry ) {
   std::vector<std::string > str(entry.size());
   for ( unsigned int j = 0 ; j<entry.size() ; j++ ) {
      std::stringstream ss;
      ss << entry[j];
      str[j]=ss.str();
   }
   AppendToTable(key,entry);
}

#endif

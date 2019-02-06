 #include "datareader.hh"

#include <stdio.h>
#include <string.h>
#include <map>
#include <string>
#include <ctype.h>
 
 
 int data_reader :: open_file(string fn, const char parser)  {
  if ( ! (file = fopen(fn.c_str(), "r")) )                return 0;
  int ret = open_file(file, parser);
  fclose(file);  
return ret;}

//-----------------------------------------------------------------------------------------
const int max_string = 4096;



int data_reader :: open_file(FILE *f, const char parser)  {
  data_reader::parser = parser;
  file = f;

  int next_str = 1;
  for (char str_main [max_string+1]; fgets(str_main, max_string, file)&&next_str;)
  {
    char  header[max_string+1], data[max_string+1];
    memset(header, 0,  max_string+1);
    memset(data,   0,  max_string+1);
    
    char *s = str_main, *h = header, *d = data;

    int dummy = 0;
    
    for (int state = 1; state;)
    {
      switch (state) {
      case 1 :
        if (*s == '#')   
        { state = 0; next_str = 0; dummy = 1;     break;}
        if (*s == parser){ s++; state = 2;        break;}        
        if (*s == ' ' || *s == '\t' || *s == '\b')   { s++; h = header;       break;}
        if (*s == 0)     { state = 0; dummy = 1;  break;}
        if (1)           { *(h++) = *(s++); *h=0; break;}
      case 2 :
        if (*s == ' ')   { s++;                   break;}
        if (1)           { state = 3;             break;}
      case 3 :
        if (*s == ' ' || *s==0 || *s=='\n' || *s == 10 || *s == 13)
        { state = 0;                              break;}
        if (1)           { *(d++) = *(s++);       break;}
      }
    }
    if (dummy) continue; 

    string s_header(header), s_data(data);
    strings.insert( std::pair<string, string>(s_header, s_data));
  }
return 1;}





//---------------------------------------------------------------------------------

int data_reader::get_string(const string &header, string &result)
{
  map <string, string>::iterator p;  
  p = strings.find(header);

  if (p != strings.end())
  {
    string s_data = p->second;
    char *str_data = new char [strlen(s_data.c_str())+1];
    if (str_data == 0) fprintf(stderr, "data_reader2.h : allocation error #3\n");    
    sscanf(s_data.c_str(), "%s", str_data) ==1? 0: str_data[0]=0;

    result = str_data;    
    delete [] str_data;    
    return 1;
  }
return 0;}

//-----------------------------------------------------------------------------------------


template <class type> type max_num(const type a, const type b) {  return a>b? a:b; }
template <class type> type min_num(const type a, const type b) {  return a<b? a:b; }

template <class type> type max_num(const type a, const type b, const type c) {  return max_num(a, max_num(b, c) );  }
template <class type> type min_num(const type a, const type b, const type c) {  return min_num(a, min_num(b, c) );  }

template <class type> type min_max_num(const type x, const type a, const type b) {  return max_num( min_num(x, b), a);  }


//==========-===============================-====================================-=========

int data_reader::get_double(const string &header, double&x){
  string data;
  if (!get_string(header, data))                  return 0;
  double tx;
  if (1 == sscanf(data.c_str(), "%lf", &tx))  {  x = tx; return 1;  }
return 0;}

//-----------------------------------------------------------------------------------------

int data_reader::get_int(const string &header, int&x){
  double t;
  if (get_double(header, t)) { x = int(t); return 1;}
return 0;}
//-----------------------------------------------------------------------------------------

int     data_reader::get_string  
(
  const string &header,
  const int     pos,
        string &result
)
{
  vector <string>  res;
  if (0 == get_list(header, res))  return 0;
  
  const int p = min_num(pos, (int)res.size()-1);
  result = res[p];
return 1;}

//-----------------------------------------------------------------------------------------

int     data_reader::get_double
(
  const string &header, 
  const int pos, 
  double &x)
{
  string data;
  if (!get_string(header, pos, data))                  return 0;
  double tx;
  if (1 == sscanf(data.c_str(), "%lf", &tx))  {  x = tx; return 1;  }
return 0;}

//-----------------------------------------------------------------------------------------

int data_reader::get_int
(
  const string &header, 
  const int pos, 
  int&x
)
{
  double t;
  if (get_double(header, pos, t)) { x = int(t); return 1;}
return 0;}      

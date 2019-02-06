#ifndef __data_reader_h__
#define __data_reader_h__

#include <map>
#include <string>
#include <vector>
#include <sstream>

  using std::map;
  using std::string;
  using std::vector;
  using std::stringstream;
  using std::istringstream;
//-----------------------------------------------------------------------------------------

  class data_reader  
  {
    FILE *file;
    map<string, string> strings;
    char parser, list_delimeter;

  public:
    data_reader ()  : list_delimeter(',')     {}
    data_reader (string fn,  const char parser='=')  : list_delimeter(',') 
      {open_file(fn,parser);}
    data_reader (FILE *f,       const char parser='=')  : list_delimeter(',') 
      {open_file(f, parser);}
    int     open_file   (FILE *f,       const char parser= '=');
    int     open_file   (string fn,     const char parser= '=');
    int     get_string  (const string &header, string &result);
    int     get_double  (const string &header, double &x);
    int     get_int     (const string &header, int &x);
    
    int     get_string  (const string &header, const int pos, string &result);
    int     get_double  (const string &header, const int pos, double &x);
    int     get_int     (const string &header, const int pos, int &x);
    

    template <class type>
    int     get_list    (const string &header, vector<type> &res)
    {
      string data;
      if (!get_string(header, data))                  return 0;

      for (int i = 0; i < (int)data.size(); i++)
      {
        data[i] = data[i] == list_delimeter? ' ' : data[i];
      }
      
      res.clear();    

      istringstream list_stream(data);

      for (type v; list_stream.good(); res.push_back(v))  list_stream >> v;

      // Ошибка чтения списка
      if(!list_stream.eof())   {  res.clear();   return 0;  }  
    return (int)res.size();}
};
//-----------------------------------------------------------------------------------------




#endif				 
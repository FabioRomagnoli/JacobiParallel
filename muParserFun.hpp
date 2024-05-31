#include <muParser.h>

#include <memory>
#include <string>

/*
functor class that relies on muparser to take as input a string and parses it.
The operator () making this a functor and making it possible to pass it to a 
funtion wrapper. The actualy evaluation is done through muparser.
*/

class MuParserFun
{
public:
  MuParserFun(const MuParserFun &m)
    : m_parser(m.m_parser)
  {
  //adress is passed to allow later evaulatiom
   m_parser.DefineVar("x", &values[0]);
   m_parser.DefineVar("y", &values[1]);
  };

  MuParserFun(const std::string &s)
  {
    try
      {
        //adress is passed to allow later evaulation
        m_parser.DefineVar("x", &values[0]);
        m_parser.DefineVar("y", &values[1]);

        m_parser.SetExpr(s);
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double
  operator()(const double &x, const double &y)
  {
    double sol = 0;

    //stored in array t hat muparser has access to
    values[0] = x;
    values[1] = y;

    try
      {
        //evulation done by muparser
        sol = m_parser.Eval();
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
    return sol;
  };

private:
  double  values[2];
  mu::Parser m_parser;
};

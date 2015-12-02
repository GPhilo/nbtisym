%{ /* -*- C++ -*- */
# include <cerrno>
# include <climits>
# include <cstdlib>
# include <string>
# include "driver.hh"
# include "parser.hh"

// Work around an incompatibility in flex (at least versions
// 2.5.31 through 2.5.33): it generates code that does
// not conform to C89.  See Debian bug 333231
// <http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=333231>.
# undef yywrap
# define yywrap() 1

// The location of the current token.
static yy::location loc;
%}
%option noyywrap nounput batch debug noinput
id    [a-zA-Z][a-zA-Z_0-9]*
float [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?
blank [ \t]

%{
  // Code run each time a pattern is matched.
  # define YY_USER_ACTION  loc.columns (yyleng);
%}

%%

%{
  // Code run each time yylex is called.
  loc.step ();
%}

{blank}+   loc.step ();
[\n]+      loc.lines (yyleng); loc.step ();
"#"[^\n]*  loc.step ();
":"		{ return yy::calcxx_parser::make_COLON(loc); }
"ID"		{ return yy::calcxx_parser::make_ID_FIELD(loc); }
"VDD"		{ return yy::calcxx_parser::make_VDD_FIELD(loc); }
"FREQUENCY"	{ return yy::calcxx_parser::make_FREQ_FIELD(loc); }
"DUTY-CYCLE" 	{ return yy::calcxx_parser::make_DC_FIELD(loc); }
{id}		{ return yy::calcxx_parser::make_IDENTIFIER(yytext, loc); }
{float}		{
		  float v = atof(yytext);
		  if( v == 0 )
			driver.error (loc, "float is invalid");
		  return yy::calcxx_parser::make_NUMBER(v, loc);
		}
.          	{ driver.error (loc, "invalid character"); }
<<EOF>>    	{ return yy::calcxx_parser::make_END(loc); }
%%

void
calcxx_driver::scan_begin ()
{
  yy_flex_debug = trace_scanning;
  if (file.empty () || file == "-")
    yyin = stdin;
  else if (!(yyin = fopen (file.c_str (), "r")))
    {
      error ("cannot open " + file + ": " + strerror(errno));
      exit (EXIT_FAILURE);
    }
}



void
calcxx_driver::scan_end ()
{
  fclose (yyin);
}


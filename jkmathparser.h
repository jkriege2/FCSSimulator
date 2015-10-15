/*
    Copyright (c) 2008-2015 Jan W. Krieger (<jan@jkrieger.de>, <j.krieger@dkfz.de>), German Cancer Research Center (DKFZ) & IWR, University of Heidelberg

    last modification: $LastChangedDate: 2015-09-30 15:27:19 +0200 (Mi, 30 Sep 2015) $  (revision $Rev: 4040 $)

    This software is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License (LGPL) as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License (LGPL) for more details.

    You should have received a copy of the GNU Lesser General Public License (LGPL)
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



/*! \defgroup jkmp parser for mathematical expressions
    \ingroup tools_math

 In this group there are classes that form a parser and evaluator for mathematical expressions.
 In the context of the sequencer program this is a tool class that can be used by the classes
 in the project. E.g. by jkiniparser which allows for mathematical expressions in the configuration
 files.

 */

/** \file jkmathparser.h
 *  \ingroup tools_math
 */

#include <iostream>
#include <cmath>
#include <string>
#include <cstdio>
#include <vector>
#include <map>
#include <sstream>
#include <exception>
#include <ctype.h>
#include <list>
#include <utility>
#include "tools.h"



#ifndef JKMATHPARSER_H
#define JKMATHPARSER_H



/**
 * \defgroup jkmpmain main function parser class
 *  \ingroup tools_math
 */
/*@{*/


/*! \brief A simple function parser to parse (build memory tree representation) and
           evaluate simple mathematical expressions
    \ingroup tools_math
 This class implements a simple function parser which can parse
 mathematical expressions like <code> z=a*3+2.34^2*sin(pi*sqrt(x))</code> .
 More than one expression can be separated by semicolon ';'. The result of
 a parse operation will be a memory structure (tree) representing the given
 expression.
 The parser can cope with constants and (user defined) variables and supports
 the data types number (double precision floating point), boolean (true/false)
 and string. It already contains a lot of fundamental mathematical functions
 (i.e. nearly all the functions from C StdLib).

 \section jkmp_constantsAndVars constants and variables:
 This class provides management utilities for constants and variables. Variables
 can also be defined as external variables, when a pointer to corresponding
 memory is provided by the calling program. The parser supports a variable
 assign operation \code a=<expression>\endcode  which allows to define new
 variables during evaluation. There are some mathematical standard constants
 registered by calling jkMathParser::addStandardVariables():
   - \c pi = \f$ \pi \f$
   - \c e = \f$ \exp(1) \f$
   - \c sqrt2 = \f$ \sqrt{2} \f$
   - \c version = the parser version
   - \c log2e = \f$ \log_2(e) \f$
   - \c log10e = \f$ \log_{10}(e) \f$
   - \c ln2 = \f$ \ln(2)=\log_{e}(2) \f$
   - \c ln10 = \f$ \ln(10)=\log_{e}(10) \f$
   - \c h = \f$ 6.6260689633\cdot 10^{-34}\;\mathrm{J\cdot s} \f$ (planck constant)
   - \c hbar = \f$ \hbar=\frac{h}{2\pi}= 1.05457162853\cdot 10^{-34}\;\mathrm{J\cdot s} \f$ (planck constant)
   - \c epsilon0 = \f$ \epsilon_0= 8.854187817\cdot 10^{-12}\;\mathrm{\frac{F}{m}} \f$ (electric constant)
   - \c mu0 = \f$ \mu_0=2\pi\cdot 10^{-7}= 12.566370614\cdot 10^{-7}\;\mathrm{\frac{N}{A^2}} \f$ (magnetic constant)
   - \c c = \f$ 299792458\;\mathrm{\frac{m}{s}} \f$ (speed of light in vacuum)
   - \c ce = \f$ 1.60217648740\cdot 10^{-19}\;\mathrm{C}\f$ (elementary charge)
   - \c muB = \f$ \mu_B= 927.40091523\cdot 10^{-26}\;\mathrm{\frac{J}{T}} \f$ (Bohr magneton)
   - \c muB_eV = \f$ \mu_B=5.788381755579\cdot 10^{-5}\;\mathrm{\frac{eV}{T}}  \f$ (Bohr magneton)
   - \c muN = \f$ \mu_N=5.0507832413\cdot 10^{-27}\;\mathrm{\frac{J}{T}}  \f$ (nuclear magneton)
   - \c muN_eV = \f$ \mu_N=3.152451232645\cdot 10^{-8}\;\mathrm{\frac{eV}{T}}  \f$ (nuclear magneton)
   - \c me = \f$ m_e= 9.1093821545\cdot 10^{-31}\;\mathrm{kg} \f$ (mass of electron)
   - \c mp = \f$ m_p= 1.67262163783\cdot 10^{-27}\;\mathrm{kg} \f$ (mass of proton)
   - \c mn = \f$ m_n= 1.67492721184\cdot 10^{-27}\;\mathrm{kg} \f$ (mass of neutron)
   - \c NA = \f$ N_A= 6.0221417930\cdot 10^{23} \f$ (Avogadro constant = particles in 1 mol)
   - \c kB = \f$ k_B=1.380650424\cdot 10^{-23}\;\mathrm{\frac{J}{K}}  \f$ (Boltzman constant)
   - \c kB_eV = \f$ k_B=8.61734315\cdot 10^{-5}\;\mathrm{\frac{eV}{K}}  \f$ (Boltzman constant)
.
 You can add user-defined contants by calling jkMathParser::addVariableDouble()
 jkMathParser::addVariableBoolean() or jkMathParser::addVariableString()


 \section jkmp_functions functions:
 this class provides a wide range of (mathematical) functions:
   - sinc, gauss, slit, theta, tanc, sigmoid
   - asin, acos, atan, atan2,
   - sin, cos, tan,
   - sinh, cosh, tanh
   - log, log2, log10,
   - exp, sqr, sqrt
   - abs, sign
   - if
   - erf, erfc, lgamma, tgamma, j0, j1, jn, y0, y1, yn
   - rand, srand
   - ceil, floor, trunc, round,
   - fmod, min, max
   - floattostr, booltostr
 .

 these functions are registere by calling jkMathParser::addStandardFunctions().
 you can add new functions by calling jkMathParser::addFunction();

 \section jkmp_resultsofparsing result of parsing and evaluation:
 The result of calling jkMathParser::parse()
 will be a tree-like structure in memory. The parse() function will return
 a pointer to the root node of this structure. All nodes inherit from
 jkmpNode class. To evaluate such a structure simply call jkmpNode::evaluate()
 of the root node. This will then return a jkmpResult structure which contains
 the result. This scheme allows for once parsing and multiply evaluating an expression.
 So if you want to define a function by an expression you can provide an
 external variable x as the argument and then evaluate the function for
 each x.

 \section jkmp_ebnf EBNF definition of the parsed expressions

<pre> logical_expression ->  logical_term
                      | logical_expression <b>or</b> logical_term
                      | logical_expression <b>||</b> logical_term</pre>
<pre> logical_term    ->  comp_expression
                      | logical_term <b>and</b> comp_expression
                      | logical_term <b>&amp;&amp;</b> comp_expression</pre>
<pre> comp_expression  ->  math_expression
                      | expression <b>==</b> math_expression
                      | expression <b>!=</b> math_expression
                      | expression <b>&gt;=</b> math_expression
                      | expression <b>&lt;=</b> math_expression
                      | expression <b>&gt;</b> math_expression
                      | expression <b>&lt;</b> math_expression</pre>
<pre>  math_expression ->  term
                      | math_expression <b>+</b> math_term
                      | math_expression <b>-</b> math_term</pre>
<pre>  math_term       ->  primary
                      | term <b>*</b> primary
                      | term <b>/</b> primary
                      | term ( <b>%</b> | <b>mod</b> ) primary</pre>
<pre>  primary         ->  <b>true</b> | <b>false</b>
                      | string_constant
                      | NUMBER
                      | NAME
                      | NAME <b>=</b> logical_expression
                      | <b>+</b> primary | <b>-</b> primary | <b>!</b> primary | <b>not</b> primary
                      | <b>(</b> logical_expression <b>)</b>
                      | NAME<b>(</b> parameter_list <b>)</b>
                      | primary <b>^</b> primary</pre>
<pre>  string_constant -> <b>&quot;</b> STRING <b>&quot;</b> | <b>&apos;</b> STRING <b>&apos;</b></pre>
<pre>  parameter_list  ->  \f$ \lambda \f$ | logical_expression | logical_expression <b>,</b> parameter_list</pre>




  \section jkmp_example Simple Example of Usage
     \code
        try {
            jkMathParser mp; // instanciate
            jkmpNode* n;
            jkmpResult r;
            // parse some numeric expression
            n=mp.parse("pi^2+4*sin(65*pi/exp(3.45))");
            r=n->evaluate();
            cout<<r.num<<endl;
            //delete n;

            // parse some boolean expression
            n=mp.parse("true==false");
            r=n->evaluate();
            if (r.type==jkmpBool) {
                if (r.boolean) cout<<"true";
                else cout<<"false";
            }
            if (r.type==jkmpDouble) cout<<r.num<<endl;
            if (r.type==jkmpString) cout<<r.str<<endl;
            delete n;

            // parse some string expression
            n=mp.parse("var1='false'; var1+'true'");
            r=n->evaluate();
            if (r.type==jkmpString) cout<<r.str<<endl;
            delete n;
        } catch(std::exception& E) {
            cout<<"ERROR!!!\n  "<<E.what()<<endl<<endl;
        }
     \endcode



     \section jkmp_errorhandling Error Handling

     In the above example we use error handling by use of exception (default behauviour).


     It is also possible to change the error handling from using exceptions to calling a specific
     error handling function. This can be usefull in programs that don't support exceptions.
     To do so, use this cod:
     \code
     void error(std::string message) {
       cout<<"error: "+message;
       system("PAUSE");
       abort();
     }

     int main() {
       jkMathParser mp;
       mp.set_exception_function(error);  // make error ahndler known
       ...
     }
     \endcode
 */
class jkMathParser
{
    protected:
        void* data;

        /** \brief the possible tokens that can be recognized by the tokenizer in jkMathParser::getToken() */
        enum jkmpTokenType {
            END,                /*!< \brief end token */
            PRINT,              /*!< \brief a semicolon ';' */
            PARAMETER_DIV,      /*!< \brief a comma ',' between two function parameters */
            STRING_DELIM,       /*!< \brief a string delimiter ' or " */
            NAME,               /*!< \brief a name (consisting of characters) of a variable or function */
            NUMBER,             /*!< \brief a number in scientific notation */
            PLUS,               /*!< \brief a plus operator '+' */
            MINUS,              /*!< \brief a minus operator '-' */
            MUL,                /*!< \brief a multiplication operator '*' */
            DIV,                /*!< \brief a division operator '/' */
            MODULO,             /*!< \brief a modulo operator '%' */
            ASSIGN,             /*!< \brief a variable assignment = */
            LBRACKET,           /*!< \brief left brackets '(' */
            RBRACKET,           /*!< \brief right brackets ')' */
            POWER,              /*!< \brief a power operator '^' */
            FACTORIAL_LOGIC_NOT, /*!< \brief a factorial operator or a logical NOT '!' */
            LOGIC_NOT,          /*!< \brief a logical NOT '!' / 'not' */
            LOGIC_AND,          /*!< \brief a logical AND operator '&&' / 'and' */
            LOGIC_OR,           /*!< \brief a logical OR operator '||' / 'or' */
            LOGIC_XOR,          /*!< \brief a logical XOR operator 'xor' */
            LOGIC_NOR,          /*!< \brief a logical NOR operator 'nor' */
            LOGIC_NAND,         /*!< \brief a logical NAND operator 'nand' */
            LOGIC_TRUE,         /*!< \brief 'true' */
            LOGIC_FALSE,        /*!< \brief  'false' */
            COMP_EQUALT,        /*!< \brief equals operation '==' */
            COMP_UNEQUAL,       /*!< \brief unequal operation '!=' */
            COMP_GREATER,       /*!< \brief greater than operation '>' */
            COMP_SMALLER,       /*!< \brief smaller than operation '<' */
            COMP_GEQUAL,        /*!< \brief greater than or equal operation '>=' */
            COMP_SEQUAL,        /*!< \brief smaller than or equal operation '<=' */
        };



        /** \brief internal names for logic operations */
        enum {
          jkmpLOPand='a',
          jkmpLOPor='o',
          jkmpLOPxor='x',
          jkmpLOPnor='n',
          jkmpLOPnand='A',
          jkmpLOPnot='n'
        };


        /** \brief jkmpCOMPdefs internal names for compare operations */
        enum {
            jkmpCOMPequal='=',
            jkmpCOMPnequal='!',
            jkmpCOMPlesser='<',
            jkmpCOMPgreater='>',
            jkmpCOMPlesserequal='a',
            jkmpCOMPgreaterequal='b'
        };
        /*@}*/

    public:
        /**
         * \defgroup jkmpultil utilities for jkMathParser function parser class
         * \ingroup jkmp
         */
        /*@{*/

        /** possible result types
         *    - \c jkmpDouble: a floating-point number with double precision. This is
         *                     also used to deal with integers
         *    - \c jkmpString: a string of characters
         *    - \c jkmpBool:   a boolean value true|false
         *  .
         */
        enum jkmpResultType {jkmpDouble,  /*!< \brief a floating-point number with double precision. This is also used to deal with integers */
                             jkmpString,  /*!< \brief a string of characters */
                             jkmpBool};   /*!< \brief a boolean value true|false */



        /** \brief result of any expression  */
        struct jkmpResult {
          jkmpResult();

          bool isValid;
          jkmpResultType type;   /*!< \brief type of the result */
          std::string str;       /*!< \brief contains result if \c type==jkmpString */
          double num;            /*!< \brief contains result if \c type==jkmpDouble */
          bool boolean;          /*!< \brief contains result if \c type==jkmpBool */

          /** \brief convert the value this struct representens into a std::string */
          inline std::string toString() {
              switch(type) {
                  case jkmpDouble: return floattostr(num);
                  case jkmpString: return str;
                  case jkmpBool: return booltostr(boolean);
              }
              return "";
          };

          /** \brief convert the value this struct representens into a std::string and adds the name of the datatype in \c [...] */
          inline std::string toTypeString() {
              switch(type) {
                  case jkmpDouble: return floattostr(num)+" [number]";
                  case jkmpString: return str+" [string]";
                  case jkmpBool: return booltostr(boolean)+" [bool]";
              }
              return "";
          }
        } ;


        /** \brief This struct is for managing variables. Unlike jkmpResult this struct
          * only contains pointers to the data
          */
        struct jkmpVariable {
          jkmpVariable();
          jkmpResultType type;     /*!< \brief type of the variable */
          bool internal;           /*!< \brief this is an internal variable */
          std::string *str;        /*!< \brief this points to the variable data if \c type==jkmpString */
          double *num;             /*!< \brief this points to the variable data if \c type==jkmpDouble */
          bool *boolean;           /*!< \brief this points to the variable data if \c type==jkmpBool */
        };

        /** \brief This struct is for managing temporary variables. It is generally like jkmpVariable.
          */
        struct jkmpTempVariable {
          std::string name;       /*!< \brief name of the variable */
          jkmpResultType type;    /*!< \brief type of the variable */
          bool internal;          /*!< \brief this is an internal variable */
          std::string *str;       /*!< \brief this points to the variable data if \c type==jkmpString */
          double *num;            /*!< \brief this points to the variable data if \c type==jkmpDouble */
          bool *boolean;          /*!< \brief this points to the variable data if \c type==jkmpBool */
        };



        /** \brief This is a function prototype for adding new mathematical functions
         *         to the parser
         *
         * If you want to add more math functions (like sin, cos , abs ...) to the
         * parser, you will have to implement it with this prototype and then register
         * it with jkMathParser::addFunction(). The first parameter points to an array
         * containing the input parameters while the second one specifies the number
         * of supplied parameters. The result has to be of type jkmpResult.
         *
         * All error handling has to be done inside the function definition. Here is a
         * simple example:
         * \code
         * jkmpResult Abs(jkmpResult* params, unsigned char n){
         *   jkmpResult r;
         *   r.type=jkmpDouble;
         *   if (n!=1) jkmpError("abs accepts 1 argument");
         *   if (params[0].type!=jkmpDouble) jkmpError("abs needs double argument");
         *   r.num=fabs(params[0].num);
         *   return r;
         * }
         * \endcode
         */
        typedef jkmpResult (*jkmpEvaluateFunc)(jkmpResult*, unsigned char, jkMathParser*);


        /** \brief description of a user registered function */
        struct jkmpFunctionDescriptor {
          jkmpEvaluateFunc function;    /*!< \brief a pointer to the function implementation */
          std::string name;             /*!< \brief name of the function */
        };


        /*@}*/


        /**
         * \defgroup jkmpNodes memory representation of expressions
         * \ingroup jkmp
         */
        /*@{*/

        /**
         * \brief This class is the abstract base class for nodes.
         * All allowed node types must inherit from jkmpNode
         */
        class jkmpNode {
          protected:
            jkMathParser* parser;  /*!< \brief points to the parser object that is used to evaluate this node */
            jkmpNode* parent;      /*!< \brief points to the parent node */
          public:
            /** \brief virtual class destructor */
            virtual ~jkmpNode() {};

            /** \brief evaluate this node */
            virtual jkmpResult evaluate()=0;

            /** \brief return a pointer to the jkMathParser  */
            inline jkMathParser* getParser(){ return parser; };

            /** \brief set the jkMathParser  */
            inline void setParser(jkMathParser* mp){ parser=mp; };

            /** \brief returns a pointer to the parent node */
            inline jkmpNode* getParent(){ return parent; };

            /** \brief sets the parent node  */
            inline void setParent(jkmpNode* par) { parent=par; };
        };


        /**
         * \brief This class represents a binary arithmetic operation:
         *  add (+), subtract (-), multiply (*), divide (/), a to the power of b (a^b)
         */
        class jkmpBinaryArithmeticNode: public jkmpNode {
          private:
            jkmpNode* left, *right;
            char operation;
          public:
            /** \brief constructor for a jkmpBinaryArithmeticNode
             *  \param op the operation to be performed: add (+), subtract (-), multiply (*), divide (/), a to the power of b (a^b)
             *  \param l left child node/operand
             *  \param r right child node/operand
             *  \param p a pointer to a jkMathParser object
             *  \param par a pointer to the parent node
             */
            jkmpBinaryArithmeticNode(char op, jkmpNode* l, jkmpNode* r, jkMathParser* p, jkmpNode* par);

            /** \brief standard destructor, also destroy the children (recursively)  */
            ~jkmpBinaryArithmeticNode() { delete left; delete right;};

            /** \brief evaluate this node */
            virtual jkmpResult evaluate();
        };

        /**
         * \brief This class represents a binary boolean operation: and, or, xor, nor, nand
         */
        class jkmpBinaryBoolNode: public jkmpNode {
          private:
            jkmpNode* left, *right;
            char operation;
          public:
            /** \brief constructor for a jkmpBinaryBoolNode
             *  \param op the operation to be performed: (a)nd, (o)r, (x)or, (n)or, nand (A)
             *  \param l left child node/operand
             *  \param r right child node/operand
             *  \param p a pointer to a jkMathParser object
             *  \param par a pointer to the parent node
             */
            jkmpBinaryBoolNode(char op, jkmpNode* l, jkmpNode* r, jkMathParser* p, jkmpNode* par);

            /** \brief standard destructor, also destroy the children (recursively)  */
            ~jkmpBinaryBoolNode() { delete left; delete right;};

            /** \brief evaluate this node */
            virtual jkmpResult evaluate();
        };

        /**
         * \brief This class represents a binary compare operation: !=, ==, >=, <=, >, <
         */
        class jkmpCompareNode: public jkmpNode {
          private:
            jkmpNode* left, *right;
            char operation;
          public:
            /** \brief constructor for a jkmpCompareNode
             *  \param op the operation to be performed: != (!), == (=), >= (b), <= (a), (>), (<)
             *  \param l left child node/operand
             *  \param r right child node/operand
             *  \param p a pointer to a jkMathParser object
             *  \param par a pointer to the parent node
             */
            jkmpCompareNode(char op, jkmpNode* l, jkmpNode* r, jkMathParser* p, jkmpNode* par);

            /** \brief standard destructor, also destroy the children (recursively)  */
            ~jkmpCompareNode () { delete left; delete right;};

            /** \brief evaluate this node */
            virtual jkmpResult evaluate();
        };

        /**
         * \brief This class represents a unary operations: ! (bool negation), - (arithmetic negation)
         */
        class jkmpUnaryNode: public jkmpNode {
          private:
            jkmpNode* child;
            char operation;
          public:
            /** \brief constructor for a jkmpUnaryNode
             *  \param op the operation to be performed: (!), (-)
             *  \param c child node/operand
             *  \param p a pointer to a jkMathParser object
             *  \param par a pointer to the parent node
             */
            jkmpUnaryNode(char op, jkmpNode* c, jkMathParser* p, jkmpNode* par);

            /** \brief standard destructor, also destroy the children (recursively)  */
            ~jkmpUnaryNode() {delete child;};

            /** \brief evaluate this node */
            virtual jkmpResult evaluate();
        };

        /**
         * \brief This class represents a variable assignment (a = expression)
         */
        class jkmpVariableAssignNode: public jkmpNode {
          private:
            jkmpNode* child;
            std::string variable;
            //char operation;
          public:
            /** \brief standard destructor, also destroy the children (recursively)  */
            ~jkmpVariableAssignNode() {delete child;};

            /** \brief constructor for a jkmpVariableAssignNode
             *  \param var name of the variable to assign to
             *  \param c child node/right-hand-side expression
             *  \param p a pointer to a jkMathParser object
             *  \param par a pointer to the parent node
             */
            jkmpVariableAssignNode(std::string var, jkmpNode* c, jkMathParser* p, jkmpNode* par);

            /** \brief evaluate this node */
            virtual jkmpResult evaluate();
        };

        /**
         * \brief This class represents a number, a string contant or a boolean contant (true/false)
         */
        class jkmpConstantNode: public jkmpNode {
          private:
            jkmpResult data;
          public:
            /** \brief constructor for a jkmpConstantNode
             *  \param d the value of the constant
             *  \param p a pointer to a jkMathParser object
             *  \param par a pointer to the parent node
             */
            jkmpConstantNode(jkmpResult d, jkMathParser* p, jkmpNode* par) { data=d; setParser(p); setParent(par); };

            /** \brief evaluate this node */
            virtual jkmpResult evaluate() { return data; };
        };

        /**
         * \brief This class represents a variable.
         */
        class jkmpVariableNode: public jkmpNode {
          private:
            std::string var;
          public:
            /** \brief constructor for a jkmpVariableNode
             *  \param name name of the variable
             *  \param p a pointer to a jkMathParser object
             *  \param par a pointer to the parent node
             */
            jkmpVariableNode(std::string name, jkMathParser* p, jkmpNode* par);

            /** \brief evaluate this node */
            virtual jkmpResult evaluate();
        };

        /**
         * \brief This class represents an arbitrary function.
         *
         * When initialized this class will get the function description that is
         * linked to the supplied function name from jkMathParser object. This
         * information is saved locally and won't be changed when evaluating!
         *
         * Functions may have 8 parameters at the most.
         */
        class jkmpFunctionNode: public jkmpNode {
          private:
            std::string fun;
            jkmpNode** child;
            unsigned char n;
            jkmpEvaluateFunc function;
          public:
            /** \brief constructor for a jkmpFunctionNode
             *  \param name name of the function
             *  \param c a pointer to an array of jkmpNode objects that represent the parameter expressions
             *  \param num number of children in c
             *  \param p a pointer to a jkMathParser object
             *  \param par a pointer to the parent node
             */
            jkmpFunctionNode(std::string name, jkmpNode** c, unsigned char num, jkMathParser* p, jkmpNode* par);

            /** \brief standard destructor, also destroy the children (recursively) */
            ~jkmpFunctionNode();

            /** \brief evaluate this node */
            virtual jkmpResult evaluate();
        };

        /**
         * \brief This class represents a list of jkmpNode.
         *
         * when evaluating the result will be the result of the last node in the list.
         */
        class jkmpNodeList: public jkmpNode {
          private:
            std::vector<jkmpNode*> list;
          public:
            /** \brief constructor for a jkmpNodeList
             *  \param p a pointer to a jkMathParser object
             */
            jkmpNodeList(jkMathParser* p) { setParser(p); setParent(NULL); };

            /** \brief standard destructor, also destroy the children (recursively) */
            ~jkmpNodeList();

            /** \brief add a jkmpNode n to the list */
            void add(jkmpNode* n);

            /** \brief evaluate the node */
            virtual jkmpResult evaluate();

            /** \brief get the number of nodes in the list */
            int getCount() {return list.size();};
        };

        /*@}*/


    public:
        /**
         * \defgroup jkmpErrorhandling error handling
         * \ingroup jkmp
         */
        /*@{*/




        /** \brief error handling: exceptions of the type of this class will be thrown if an error occurs
         *
         * \attention If you do not want to use the exception handling which throws
         * jkmpException exceptions, but want to write your own error handling, you should write your own
         * error handler and assign it (function pointer) to the global variable jkmathparser_exception_function.
         * If this is not NULL this function will be called instead of throwing an exception.
         */
        class jkmpException : public std::exception {
          private:
               /** \brief the error message */
               std::string errormessage;
            public:
                /** \brief class constructors */
                jkmpException() throw() {
                    errormessage="unknown error";
                };

                /** \brief constructor with supplied error message */
                jkmpException(std::string msg) throw() {
                    errormessage=msg;
                };

                /** \brief class destructors */
                ~jkmpException() throw() {   };

                /** \brief  returns the assigned errormessage */
                inline std::string getMessage() const {
                    return errormessage;
                };

                /** \brief returns the error description as C string */
                virtual const char* what() const throw() { return getMessage().c_str(); };
        };

        /** \brief type for a custom error handler. This an alternative error handling ... may be used together with Matlab in a MEX file! */
        typedef void (*jkmpexceptionf)(std::string);


        /** \brief macro that throws an exception or calls an error handler */
        inline void jkmpError(std::string st) {
            if (jkmathparser_exception_function!=NULL) {
                jkmathparser_exception_function(st);
            } else {
                throw jkmpException(st);
            }
        };

    private:
        /** \brief if this is NULL then an exception may be thrown otherwise this should point to an error handler that will be called. */
        jkmpexceptionf jkmathparser_exception_function;

    public:
        /** \brief activate error handling by use of an exception function */
        inline void set_exception_function(jkmpexceptionf exception_function) {
            jkmathparser_exception_function=exception_function;
        };

        /** \brief deactivate error handling by use of an exception function */
        inline void reset_exception_function() {
            jkmathparser_exception_function=NULL;
        };


        /*@}*/

    protected:
        /** \brief return the given token as human-readable string */
        std::string tokentostring(jkmpTokenType token);

        /** \brief return the current token as human-readable string */
        std::string currenttokentostring();

		/** \brief Tokenizer: extract the next token from the input */
		jkmpTokenType getToken();

		/** \brief return a delimited text, i.e. extract the texte between the delimiters <code>"</code> in: of <code>"Hallo!"</code>, i.e. returns <code> Hallo!</code>
         *         This is used to parse string constants.
         *
         * This functions actually reads pascal style delimited string constants. So if you want to use the delimiter as part of the string you will have to
         * write it as doubled character. So <code>'Jan''s Test'</code> stands for <code>Jan's Test</code>.
         */
		std::string readDelim(char delimiter);

		/** \brief recognizes an compExpression while parsing. If \a get ist \c true, this function first retrieves a new token by calling getToken() */
		jkmpNode* compExpression(bool get);

		/** \brief recognizes a logicalExpression while parsing. If \a get ist \c true, this function first retrieves a new token by calling getToken() */
		jkmpNode* logicalExpression(bool get);

		/** \brief recognizes a logicalTerm while parsing. If \a get ist \c true, this function first retrieves a new token by calling getToken() */
		jkmpNode* logicalTerm(bool get);

		/** \brief recognizes a mathExpression while parsing. If \a get ist \c true, this function first retrieves a new token by calling getToken() */
		jkmpNode* mathExpression(bool get);

		/** \brief recognizes a term while parsing. If \a get ist \c true, this function first retrieves a new token by calling getToken() */
		jkmpNode* mathTerm(bool get);

		/** \brief recognizes a primary while parsing. If \a get ist \c true, this function first retrieves a new token by calling getToken() */
		jkmpNode* primary(bool get);

		/** \brief this stream is used to read in the program. An object is created and assigned
  		 * (and destroyed) by the parse()-function */
		std::istringstream* program;

        /** \brief vector containing all temporary variables */
        std::vector<jkmpTempVariable> tempvariables;

        /** \brief map to manage all currently defined variables */
        std::map<std::string, jkmpVariable> variables;

        /** \brief map to manage all currently rtegistered functions */
        std::map<std::string, jkmpFunctionDescriptor> functions;

        /** \brief the current token while parsing a string */
        jkmpTokenType CurrentToken;

        /** \brief the string value of the current token (when applicable) during the parsing step */
        std::string StringValue;

        /** \brief the string value of the current token (when applicable) during the parsing step */
        double NumberValue;

        /** \brief set the defining struct of the given variable */
        void setVariable(std::string name, jkmpResult value);
        /** \brief set the defining struct of the given variable */
        void setVariableDouble(std::string name, double value);


        /**\brief  adds a temporary variable */
        void addTempVariable(std::string name, jkmpResult value);

	protected:
		int argc;
		char **argv;

	public:
		/**\brief class constructor */
		jkMathParser();

		/**\brief class destructor */
		virtual ~jkMathParser();

        inline void* get_data() const { return data;}
        inline void set_data(void* d)  { data=d;}

		/**\brief  register a new function
		 * \param name name of the new function
		 * \param function a pointer to the implementation
		 */
		void addFunction(std::string name, jkmpEvaluateFunc function);

		/**\brief  register a new external variable of type double
		 * \param name name of the new variable
		 * \param v pointer to the variable memory
		 */
		void addVariableDouble(std::string name, double* v);

		/**\brief  register a new external variable of type string
		 * \param name name of the new variable
		 * \param v pointer to the variable memory
		 */
		void addVariableString(std::string name, std::string* v);

		/**\brief  register a new external variable of type boolean
		 * \param name name of the new variable
		 * \param v pointer to the variable memory
		 */
		void addVariableBoolean(std::string name, bool* v);


		/**\brief  register a new internal variable of type double
		 * \param name name of the new variable
		 * \param v initial value of this variable
		 */
		void addVariableDouble(std::string name, double v);

		/**\brief  register a new internal variable of type string
		 * \param name name of the new variable
		 * \param v initial value of this variable
		 */
		void addVariableString(std::string name, std::string v);

        /**\brief  register a new internal variable of type boolean
         * \param name name of the new variable
         * \param v initial value of this variable
         */
        void addVariableBoolean(std::string name, bool v);

        /**\brief  register a new internal variable of type boolean
         * \param name name of the new variable
         * \param v initial value of this variable
         */
        void addVariable(std::string name, jkmpResult result);



        /**\brief  returns the value of the given variable */
        jkmpResult getVariable(std::string name);
        /**\brief  returns the value of the given variable */
        jkmpResult getVariableOrInvalid(std::string name);

        /**\brief  returns the defining structure of the given variable */
        jkmpVariable getVariableDef(std::string name);


        /**\brief  evaluates a registered function
         * \param name name of the (registered function) to be evaluated
         * \param params array of the input parameters
         * \param n number of input parameters (<=8)
         */
        jkmpResult evaluateFunction(std::string name, jkmpResult* params, unsigned char n);

        /**\brief  returns the defining structure of the given function */
        jkmpEvaluateFunc getFunctionDef(std::string name);

        /**\brief  tests whether a temporary variable exists */
        inline bool tempvariableExists(std::string name){
          if (tempvariables.size()<=0)  return false;
          for (int i=tempvariables.size()-1; i>=0; i--) {
            if (tempvariables[i].name==name) return true;
          }
          return false;
        };

        /**\brief  tests whether a variable exists */
        inline bool variableExists(std::string name){ return tempvariableExists(name)||(variables.find(name)!=variables.end()); };

        /**\brief  tests whether a function exists */
        inline bool functionExists(std::string name){ return !(functions.find(name)==functions.end()); };

        /**\brief  deletes all defined variables. the memory of internal variables
         * will be released. the external memory will not be released.
         */
        void clearVariables();

        /**\brief  delete the specified variabale and releases its internal memory.*/
        void deleteVariable(std::string name);

        /**\brief  clears the list of internal functions*/
        inline void clearFunctions() {functions.clear();};

        /**\brief  registers standard variables*/
        void addStandardVariables();

        /**\brief  registers standard functions*/
        void addStandardFunctions();

        /**\brief  parses the given expression*/
        jkmpNode* parse(std::string prog);

        /** \brief evaluate the given expression */
        jkmpResult evaluate(std::string prog);

        /**\brief  prints a list of all registered variables */
        void printVariables();

        std::vector<std::pair<std::string, jkmpVariable> > getVariables();

		void setArgCV(int argc, char **argv);

		std::string getArgCVParam(std::string name, std::string defaultResult);
};

#endif // JKMATHPARSER_H
/*@}*/

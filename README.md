# JohansenCointegration
A C++ implementation of the Johansen Cointegration test

This is the first C/C++ implementation of Johansen's Cointegration test I'm aware of on the public Internet.  It's based on a C# implementation from a poster named <a href="https://web.archive.org/web/20160727212132/http://www.quantcode.com/userinfo.php?uid=190">Vanna</a> from quantcode.com.

This isn't the best code I've written but in my defense the algorithm is very involved and the original C# code was written by someone else.  However it does in fact work and even has all the memory leaks worked out.

Prerequisites:

This class is based on the latest GNU Scientific Library (GSL) version.  If you're on a *nix system this shouldn't be a problem.  Just download and install it the usual way for your system.  Google is your friend.

For Windows systems this is stickier.  GSL does not officially provide library binaries for Windows;  you'll either have to find a Windows binary package somewhere on the Interwebs or compile it yourself (not for the faint of heart).  I successfully created my own Windows GSL binaries using the Github repository <a href=https://github.com/BrianGladman/gsl>here</a>.  The repo provides the instructions to build GSL in Visual Studio.

Build instructions for the Johansen Cointegration project:

On Linux, create a build subdirectory and go into it.  Type "cmake .." followed by make.

On Windows, create a build subdirectory, go into it and type "cmake ..".  This creates the solution and project file for Visual Studio.  The Visual Studio project is expecting your Windows-based GSL library to be in the top level directory of this project, parallel with the source files.

The main class for doing all the Johansen work is in JohansenHelper.h/.cpp.  A command line program is contained in JohansenTest.cpp which uses the JohansenHelper.  To demonstrate the cointegration function there are 4 test cases implemented in JohansenTest.  These cases were covered in Vanna's quantcode postings on quantcode.  S/he wrote many C# apps and online papers/FAQs which I encourage you to have a look at.

In particular, I referred to this article frequently: <a href="https://web.archive.org/web/20161025043033/http://www.quantcode.com/modules/smartfaq/faq.php?faqid=103">How do I interpret Johansens' test results?</a>  This is the article from where I got the 4 test cases.

Drop a line if you find this useful or have problems compiling or using this.

# JohansenCointegration
A C++ implementation of the Johansen Cointegration test

This is the first C/C++ implementation of Johansen's Cointegration test I'm aware of on the public Internet.  It's based on a C# implementation from a poster named <a href=http://www.quantcode.com/userinfo.php?uid=190>Vanna</a> from quantcode.com.

Prerequisites:

This class is based on the GNU Scientific Library (GSL) version 1.16.  If you're on a *nix system this shouldn't be a problem.  Just download and install it the usual way for your system.  Google is your friend.

For Windows systems this is stickier.  GSL does not officially provide library binaries for Windows;  you'll either have to find a Windows binary package somewhere on the Interwebs or compile it yourself (not for the faint of heart).  <a href=http://brgladman.org/oldsite/computing/gnu_scientific_library.php>This website</a> provided all the tools and detailed instructions necessary to do this.  In fact I successfully created my own Windows GSL binaries using his instructions.

Build instructions:

On Linux, create a build subdirectory and go into it.  Type "cmake .." followed by make.

On Windows, I uploaded the .sln and .vcxproj files for a Visual Studio 2013 project.  Eventually I'll fix the CMakeLists.txt file to generate VS project files directly.  The Visual Studio project is expecting your Windows-based GSL library to be in the top level directory of this project, parallel with the source files.

The main class for doing all the Johansen work is in JohansenHelper.h/.cpp.  A command line program is contained in JohansenTest.cpp which uses the JohansenHelper.  To demonstrate the cointegration function there are 4 test cases implemented in JohansenTest.  These cases were covered in Vanna's quantcode postings on quantcode.  S/he wrote many C# apps and online papers/FAQs which I encourage you to have a look at.

In particular, I referred to this article frequently: <a href=http://www.quantcode.com/modules/smartfaq/faq.php?faqid=103>How do I interpret Johansens' test results?</a>

Drop a line if you find this useful or have problems compiling or using this.

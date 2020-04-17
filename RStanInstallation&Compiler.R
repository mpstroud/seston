

#I was able to make the compiler work on the mac follwoing intrsuctions from this reposetory:
https://github.com/stan-dev/rstan/wiki/RStan-Mojave-Mac-OS-X-Prerequisite-Installation-Instructions


#This page describes the prerequisite C++ toolchain installation procedure for RStan on the Mojave version Mac OS X. They should also work for pre-Mojave versions of Mac OS X, but in that case we instead recommend a different installer. Mojave seems to break a lot of things in different ways for different people, so if you run into trouble, ask on https://discourse.mc-stan.org/ .

#First, determine if you already have Xcode installed

#Type command + space to open Spotlight.
#In the upper right of your screen, where Spotlight opens, enter xcode into the search box.
#If Xcode does not pop up under "Applications" heading, you do not have Xcode installed.
#Open "App Store" application
#Search for 'xcode', choose Xcode, and follow its steps to install Xcode (using all default options is fine for Stan).
#Open Xcode once to accept the license agreement (don't forget this step!!)
#If Xcode does pop up under "Applications", make sure you have the latest version.
#Click command + space to open Spotlight.
#Enter enter app store into search box for Spotlight in upper right of your screen.
#Click on App Store under Applications to open.
#In App Store, click on Updates menu item (below blue downarrow).
#If Xcode has an update, click to install.
#Open Xcode once to accept the license agreement (don't forget this step!!)

#Thanks to the notes here by @coatless for the instructions below. First, run the command line in terminal:
  
  sudo installer -pkg \
/Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg \
-target /

#Then, fix ~/.R/Makevars:
  
  sudo touch ~/.R/Makevars && sudo echo 'CPPFLAGS="-isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include"' >> ~/.R/Makevars

#--------- Then start R and run install.packages('rstan').

remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")

install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

#check the C++ toolchain
pkgbuild::has_build_tools(debug = TRUE)

library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores())


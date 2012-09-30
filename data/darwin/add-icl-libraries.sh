#!/bin/bash
cp /opt/intel/composer_xe_*/compiler/lib/libiomp5.dylib Mitsuba.app/Contents/Frameworks
find Mitsuba.app/Contents/MacOS/ Mitsuba.app/plugins -type f | xargs -n 1 install_name_tool -change libiomp5.dylib @rpath/libiomp5.dylib
find Mitsuba.app/Contents/Frameworks/libmitsuba-* -type f | xargs -n 1 install_name_tool -change libiomp5.dylib @rpath/libiomp5.dylib

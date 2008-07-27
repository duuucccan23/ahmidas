INCLUDE_DIRECTORIES(BEFORE SYSTEM /home/albert/Development/Ahmidas/trunk/src)
SET(Sources_clime c-lime/lime_fseeko.cpp c-lime/lime_header.cpp c-lime/lime_reader.cpp c-lime/lime_utils.cpp c-lime/lime_writer.cpp)
SET(Sources_Reader Reader/Reader_Reader.cpp Reader/Reader_read.cpp Reader/Reader_~Reader.cpp Reader/Reader_retrieveRecord.cpp)
ADD_LIBRARY(Lime ${Sources_clime} ${Sources_Reader})

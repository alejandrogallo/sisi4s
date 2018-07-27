set(CTF "${PROJECT_SOURCE_DIR}/lib/build/${CONFIG}/ctf")
find_path (
  CTF_DIR lib/libctf.a
  HINTS ENV CTF_DIR
  PATHS ${CTF} "/usr/lib/ctf" "/usr/local/lib/ctf"
  DOC "CTF Directory"
)


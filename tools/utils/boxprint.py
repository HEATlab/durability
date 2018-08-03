##
# \file boxprint.py
#
# \brief Utility file for pretty box printing.
#
# \details
#   Prints a message like so:
# \verbatim
#-------
#Message
#-------
# \endverbatim
#

## Pretty print a message in a box.
def boxprint(message, boxtopchar='-', boxsidestring=None):

  mess_len = len(message)

  # Check that the arguments will not work.
  if len(boxtopchar) != 1 or not(type(boxtopchar) is str):
    raise ValueError(
      "Incorrect boxtopchar for boxprint")

  if boxsidestring != None:
    if not(type(boxsidestring) is str):
      raise ValueError(
        "Incorrect boxsidechar for boxprint")


  if boxsidestring != None:
    print(boxtopchar*(mess_len+2*len(boxsidestring)+2))
    print(boxsidestring+" "+message+" "+boxsidestring)
    print(boxtopchar*(mess_len+2*len(boxsidestring)+2))
  else:
    print(boxtopchar*mess_len)
    print(message)
    print(boxtopchar*mess_len)


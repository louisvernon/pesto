# this module allows shared variables to be accessed across modules.
# While I don't really like this practice in general, it's useful for the status variable


import Queue

status_queue = Queue.Queue()
status = ""
import smtplib

class Gmail(object):
	SERVER = "smtp.gmail.com"
	PORT = 587
	
	def __init__(self, user, password):
		self.user = user
		self.password = password
	
	def send(self, to, subject, msg):
		if isinstance(to, list):
			to = ",".join(to)
			
		server = smtplib.SMTP(Gmail.SERVER, Gmail.PORT)
		server.ehlo()
		server.starttls()
		server.ehlo()
		server.login(self.user, self.password)
		
		payload = ["To:{0}".format(to)]
		payload += ["From:{0}".format(self.user)]
		payload += ["Subject:{0}".format(subject)]
		payload += ["\n{0}\n\n".format(msg)]
		
		server.sendmail(self.user, to, "\n".join(payload))
		
		server.close()

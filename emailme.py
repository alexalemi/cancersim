#Python script to send simple emails
import smtplib
from email.mime.text import MIMEText
import netrc

DEFAULTFROM = "alexalemi@gmail.com"
DEFAULTTO = "alexalemi@gmail.com"
DEFAULTSUBJECT = "PYTHON NOTIFICATION"

"""
Note that by default, this looks to check your netrc credentials
to use this feature, create a .netrc file, so that only you can read and write it

    touch ~/.netrc
    chmod 600 ~/.netrc

and then add the information for the gmail smtp server, i.e.

    machine smtp.gmail.com
            login yourusername@gmail.com
            password yourpassword

This way only you will have access to this file
"""

def send_email(subject=DEFAULTSUBJECT,
    message=".",me=DEFAULTFROM,recipients=[DEFAULTTO],
    smtpserver="smtp.gmail.com",tls=True,login=None,
    password=None):
    """ Send an email using the gmail smtp servers, and netrc to hide the username
        and password information """

    if login is None or password is None:
        secrets = netrc.netrc()
        netrclogin,netrcaccount,netrcpassword = secrets.authenticators(smtpserver)
    if login is None:
        login = netrclogin
    if password is None:
        password = netrcpassword

    msg = MIMEText(message)

    msg['Subject'] = subject
    msg['From'] = me
    msg['To'] = ", ".join(recipients)

    s = smtplib.SMTP(smtpserver)
    if tls:
        s.starttls()
        s.login(login,password)

    s.sendmail(me,recipients, msg.as_string())
    s.quit()


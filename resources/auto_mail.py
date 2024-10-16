# create a scipt that sends an email to a list of recipients

# Imports:
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import pandas as pd
import smtplib
from dataclasses import dataclass


# Function to send email
@dataclass
class AutoEmail():
    '''Class to send email'''

    def __init__(self, email, password, subject, message, recipients):
        self.email = email
        self.password = password
        self.subject = subject
        self.message = message
        self.recipients = recipients
    

    def read_excel(self, file):
        pass

    def send_email(self):
        pass

    def email_body(self):
        pass

    def process_email(self):
        pass

    


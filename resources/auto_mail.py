# create a scipt that sends an email to a list of recipients

# Imports:
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from dataclasses import dataclass
import pandas as pd
import smtplib

#! Remoove before deployment 
from icecream import ic

# Function to send email
@dataclass
class AutoEmail():
    '''Class to send email'''   

    def read_excel(self, file):  
        """
        Reads an Excel file and removes rows with missing or invalid data.

        Conditions for removal:
        - 'Status' is not NaN (already contacted).
        - Both 'First Name' and 'Last Name' are NaN.
        - 'Email' is NaN.
        - 'Company' is NaN.

        Args:
            file (str): Path to the Excel file.

        Modifies:
            self.data (pandas.DataFrame): The DataFrame with invalid rows removed.
        """  
        
        self.data = pd.read_excel(file)
        dropped_rows = []

        # Handle all of the edge cases
        for index, row in self.data.iterrows():
            
            # if we have already reached out to the company, skip them
            if pd.notna(row['Status']):
                ic(f"Already reached out to: {row['Company']}")
                dropped_rows.append(index)
            
            # If we are missing a name for the contact, skip them
            if pd.isna(row['First Name']) and pd.isna(row['Last Name']):
                ic(f"Missing full name for: {row['Company']}")
                dropped_rows.append(index)
            
            # If we are missing the email, skip them
            if pd.isna(row['Email']):
                ic(f"Missing email for: {row['Company']}")
                dropped_rows.append(index)

            # If we are missing the company name, skip them
            if pd.isna(row['Company']):
                ic(f"Missing company name at Row: {index + 2}")
                dropped_rows.append(index)

        self.data.drop(dropped_rows, inplace=True)
        

    def send_email(self):
        
        self.full_msg = f'''Dear {self.first_name} {self.last_name}, \n\n{self.msg} \n\nKindly, '''

    def email_body(self):
        
        self.msg = f"""
                    I hope you're well. My name is Winston Price, Team Lead for Solaris Propulsion at Embry-Riddle Aeronautical University. Our team is developing an aerospike engine designed to optimize performance across varying altitudes. Our goal is to successfully conduct two hot fire tests, achieving 1,800 pounds of force for 10 seconds, and to aid future teams in their efforts to surpass the Kármán line.

In addition to the engine, we are developing software tools to support the design and optimization of future aerospike engines, providing valuable resources for the aerospace community.

We are currently seeking sponsorship to help 3D print our engine, a critical step that will enable us to leverage additive manufacturing for greater design efficiency and cost-effectiveness.

I would be happy to discuss potential collaboration in more detail. Please let me know if you're interested in learning more or scheduling a meeting.

Thank you for your consideration."""

    def process_email(self):
        pass

if __name__ == '__main__':
    email = AutoEmail()
    email.read_excel('/Users/winston/Desktop/School/Classes/Senior 1/Capstone/Solaris-Propulsion/resources/Test Companies.xlsx')

    


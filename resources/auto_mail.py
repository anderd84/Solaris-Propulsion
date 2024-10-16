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

    def __init__(self, file: str):

        self.from_email = 'solaris.propulsion@gmail.com'
        self.from_password = 'tjql josh ykvm ohev' # App password do not change
        self.smtp_server = 'smtp.gmail.com'
        self.smtp_port = 587


        self.read_excel()
        self.email_body()
        self.send_email()
        self.update_excel()

    
    def read_excel(self):  
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
                print(f"Already reached out to: {row['Company']}")
                dropped_rows.append(index)
            
            # If we are missing a name for the contact, skip them
            if pd.isna(row['First Name']) and pd.isna(row['Last Name']):
                print(f"Missing full name for: {row['Company']}")
                dropped_rows.append(index)
            
            # If we are missing the email, skip them
            if pd.isna(row['Email']):
                print(f"Missing email for: {row['Company']}")
                dropped_rows.append(index)

            # If we are missing the company name, skip them
            if pd.isna(row['Company']):
                print(f"Missing company name at Row: {index + 2}")
                dropped_rows.append(index)

        self.data.drop(dropped_rows, inplace=True)
        

    def send_email(self):

        """
        Sends an email to each valid contact in the filtered self.data DataFrame.
        The email body is dynamically generated based on the contact's details.
        """

        for index, row in self.data.iterrows():
            self.first_name = row['First Name']
            self.last_name = row['Last Name']
            self.company = row['Company']
            recipient_email = row['Email']
            
            # Composing the message
            self.full_msg = f'''Dear {self.first_name} {self.last_name}, \n\n{self.msg} \n\nKindly,\nWinston Price \nTeam Lead, Solaris Propulsion\nEmbry-Riddle Aeronautical University'''

            # Create the email object
            msg = MIMEMultipart()
            msg['From'] = self.from_email
            msg['To'] = recipient_email
            msg['Subject'] = 'Outreach and Sponsorship Request'

            # Attach the body text to the email
            msg.attach(MIMEText(self.full_msg, 'plain'))

            try:
                # Set up the SMTP server
                server = smtplib.SMTP(self.smtp_server, self.smtp_port)
                server.starttls()  # Secure the connection
                server.login(self.from_email, self.from_password)

                # Send the email
                server.send_message(msg)
                print(f"Email successfully sent to {recipient_email}.")

            except Exception as e:
                print(f"Failed to send email to {recipient_email}. Error: {str(e)}")

            finally:
                server.quit()

    
    def email_body(self):
        
        self.msg = f"""My name is Winston Price, Team Lead for Solaris Propulsion at Embry-Riddle Aeronautical University. Our team is developing an aerospike engine designed to optimize performance across varying altitudes. Our goal is to successfully conduct two hot fire tests, achieving 1,800 pounds of force for 10 seconds, and to aid future teams in their efforts to surpass the Kármán line.

In addition to the engine, we are developing software tools to support the design and optimization of future aerospike engines, providing valuable resources for the aerospace community.

We are currently seeking sponsorship to help 3D print our engine, a critical step that will enable us to leverage additive manufacturing for greater design efficiency and cost-effectiveness.

I would be happy to discuss potential collaboration in more detail. Please let me know if you're interested in learning more or scheduling a meeting.

Thank you for your consideration."""

    
    def update_excel(self):
        """
        Updates the 'Status' column in the Excel file to 'Sent' for successfully contacted recipients.
        This will only modify the 'Status' column and leave the rest of the file intact.
        """
        
        # Load the original Excel file
        original_data = pd.read_excel(file)

        # Iterate through the data to update the 'Status' column where the email was successfully sent
        for index, row in self.data.iterrows():
            # Find the row in the original Excel that matches the email in the current DataFrame
            if pd.notna(row['Email']) and row['Email'] in original_data['Email'].values:
                original_data.loc[original_data['Email'] == row['Email'], 'Status'] = 'Sent'
        
        # Save the updated DataFrame back to the Excel file, ensuring it doesn't overwrite other data
        original_data.to_excel(file, index=False)
        print(f"Excel file {file} updated with 'Sent' status for applicable rows.")



if __name__ == '__main__':

    file = '/Users/winston/Desktop/School/Classes/Senior 1/Capstone/Solaris-Propulsion/resources/Test Companies.xlsx'
    email = AutoEmail(file)

    


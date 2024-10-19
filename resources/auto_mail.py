# Imports:
from email.mime.application import MIMEApplication
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
        self.signature_file_path = 'resources/signature.html'

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
        

    def read_signature(self):
        """
        Reads the email signature from the HTML file.
        """
        try:
            with open(self.signature_file_path, 'r') as f:
                self.signature_html = f.read()
        except Exception as e:
            print(f"Error reading signature file: {e}")
            self.signature_html = ""

    
    def send_email(self):
        """
        Sends an email to each valid contact in the filtered self.data DataFrame.
        The email body includes a manually appended signature read from an external file.
        Attach a file to the email as well.
        """

        # File you want to attach
        attachment_file_path = 'resources/Solaris Propulsion Overview.pdf'  # Update with your actual file path

        # Read the signature from the external file
        self.read_signature()

        for index, row in self.data.iterrows():
            self.first_name = row['First Name']
            self.last_name = row['Last Name']
            self.company = row['Company']
            recipient_email = row['Email']

            # Composing the message (with the manually appended signature)
            self.full_msg = f'''<p>Dear {self.first_name} {self.last_name},</p>
                                <p>{self.msg}</p> {self.signature_html}'''

            # Create the email object
            msg = MIMEMultipart('related')
            msg['From'] = self.from_email
            msg['To'] = recipient_email
            msg['Subject'] = 'Outreach and Sponsorship Request'

            # Create the body with both plain text and HTML alternatives
            msg_alternative = MIMEMultipart('alternative')
            msg.attach(msg_alternative)

            # Add the plain text version (for email clients that don't support HTML)
            msg_text = MIMEText(f'Dear {self.first_name} {self.last_name},\n\n{self.msg}\n\nKindly,\nWinston Price', 'plain')
            msg_alternative.attach(msg_text)

            # Add the HTML version
            msg_html = MIMEText(self.full_msg, 'html')
            msg_alternative.attach(msg_html)

            # Attach the file
            try:
                with open(attachment_file_path, 'rb') as attachment:
                    part = MIMEApplication(attachment.read(), Name='Solaris Propulsion.pdf') 
                    part['Content-Disposition'] = f'attachment; filename="Solaris Propulsion.pdf"'
                    msg.attach(part)
            except Exception as e:
                print(f"Failed to attach file: {e}")

            try:
                # Set up the SMTP server
                server = smtplib.SMTP(self.smtp_server, self.smtp_port)
                server.starttls()  # Secure the connection
                server.login(self.from_email, self.from_password)

                # Send the email
                server.send_message(msg)
                print(f"Email successfully sent to {recipient_email}")

            except Exception as e:
                print(f"Failed to send email to {recipient_email}. Error: {str(e)}")

            finally:
                server.quit()


    
    def email_body(self):
        
        self.msg = f"""
<p>My name is Winston Price, and I am the Team Lead for Solaris Propulsion, a Capstone team at Embry-Riddle Aeronautical University. We are currently developing an aerospike engine designed to optimize performance across varying altitudes. Our objective is to successfully complete two hot fire tests, achieving 1,800 pounds of thrust for 10 seconds, which will push the boundaries of current propulsion technology. Additionally, we aim to set the foundation for future teams to surpass the Kármán line at 330,000 feet.</p>

<p>Alongside our engine development, we are also creating software tools to assist with the design and optimization of future aerospike engines. These tools will provide valuable resources to the broader aerospace community, enabling more efficient and informed propulsion system design.</p>

<p>To make our vision a reality, we are seeking sponsorship to support the 3D printing of our engine. Utilizing additive manufacturing will allow us to achieve greater design flexibility and cost savings, advancing both our project and the future of aerospike technology.</p>

<p>I have attached a document that provides further details about our team and our mission. I would love the opportunity to discuss how your support could contribute to our success. Please let me know if you're interested in exploring potential collaboration or scheduling a meeting.</p>

<p>Thank you for your time and consideration. I look forward to hearing from you.</p>

<p>Be a part of the Future. Be a part of Solaris Propulsion.</p>

<p>Kindly,</p>
"""



    
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
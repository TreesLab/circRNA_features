import re
import urllib.request
import os.path

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC


class RBPmapJob:
    _result_url_pat = r'http://rbpmap\.technion\.ac\.il/([0-9]+)/results\.html'
    RBPmap_url = "http://rbpmap.technion.ac.il/index.html"

    def __init__(self, driver, job_name, fa_file, high_cons=False):
        self._driver = driver
        self.job_name = job_name
        self.fa_file = os.path.abspath(fa_file)
        self.high_cons = high_cons
        
        self._job_id = None
    
        self._wait = WebDriverWait(self._driver, 10)

    def init_page(self):
        if self._driver.current_url != self.RBPmap_url:
            self._driver.get(self.RBPmap_url)

    def input_data(self):
        
        # select 'upload a file'
        self._driver.find_element_by_xpath('//input[@value="file"]').click()
        
        # Motif db selection
        # open db selection page
        self._driver.find_elements_by_xpath('//input[@name="db_selection"]')[1].click()
        self._driver.find_element_by_xpath('//a[@href="RBPmap_motif_selection.html"]').click()
        
        # switch to db selection page
        self._driver.switch_to.window(self._driver.window_handles[-1])
        
        # choose the needed db
        # self._driver.find_element_by_xpath('//input[@name="human_motifs_checkbox"]').click()
        self._wait.until(EC.element_to_be_clickable((By.XPATH, '//input[@name="human_motifs_checkbox"]'))).click()
        self._wait.until(EC.presence_of_element_located((By.XPATH, '//input[@name="drosophila_motifs_checkbox"]')))

        # submit the selection
        self._driver.execute_script('SendData();')
        
        # switch back to the main page
        self._driver.switch_to.window(self._driver.window_handles[-1])

        # set "high stringency" and apply "conservation filter"
        if self.high_cons:
            self._driver.find_element_by_id('advanced_img').click()

            stringency_select = Select(self._driver.find_element_by_xpath('//select[@name="stringency"]'))
            stringency_select.select_by_value('high')

            self._driver.find_element_by_xpath('//input[@name="is_conservation"]').click()
        
        # upload file
        self._driver.find_element_by_xpath('//input[@name="seq_file"]').send_keys(self.fa_file)

    def submit(self):
        self._driver.find_element_by_xpath('//input[@value="Submit"]').click()
        
        self._wait.until(EC.url_matches(self._result_url_pat))
        self._job_id = self._get_job_id(self._driver.current_url)
    
    @property
    def job_id(self):
        return self._job_id
    
    @classmethod
    def _get_job_id(cls, url):
        job_id = re.search(cls._result_url_pat, url)[1]
        return job_id
    
    @property
    def status(self):
        return self.get_job_status(self._job_id)
    
    @staticmethod
    def get_job_status(job_id, timeout=60):
        if not job_id:
            return 'Job not submitted yet!'
        
        url = f'http://rbpmap.technion.ac.il/{job_id}/results.html'

        html_text = urllib.request.urlopen(url, timeout=timeout).read().decode()

        job_status = re.search(r'Job status: <span [^<>]*>(\w+)?</span>', html_text)[1]

        return job_status
    
    def save_result(self, out_file):
        self._save_result(self._job_id, out_file)

    @staticmethod
    def _save_result(job_id, out_file):
        if not job_id:
            return 'Job not submitted yet! Nothing to save!'
            
        result_file_url = f'http://rbpmap.technion.ac.il/{job_id}/All_Predictions.txt'
        
        with urllib.request.urlopen(result_file_url) as f_in, \
                open(out_file, 'w') as out:
            
            out.write(f_in.read().decode())














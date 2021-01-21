import os
import json
import math
from .color_space import CVC

class Litmus():
    db = []
    cell = {}
    supernova = []
    family = {}

    @classmethod
    def initialize(cls, method):
        if method == "Json":
            
            with open("static/secret/LitmusDB 20210115.json") as f:
                dj = json.loads(f.read())
            
            for index, value in enumerate(dj):
                name = value['Name']
                hexa = value['Hexa']
                RGB = CVC.hexa_RGB(hexa)
                rgb = CVC.hexa_rgb(hexa)
                cell = CVC.RGB_Cell(RGB)
                Ypq = CVC.rgb_Ypq(rgb)
                section = CVC.Ypq_Section(Ypq)
                wheel = CVC.Section_Wheel(section)['Wheel']
                depth = CVC.Section_Wheel(section)['Depth']
                mode = CVC.Section_Wheel(section)['Mode']
                group = CVC.Section_Group(section)
                star = 6 - int(value['Class'][0:1])
                text = cls.get_text(mode)
                
                litmus = {
                    'id': int(index), 
                    'name': name, 
                    'hexa': hexa,
                    'star': star,
                    'rgb': rgb,
                    'cell': cell,
                    'section': section,
                    'group': group,
                    'wheel': wheel,
                    'depth': depth,
                    'mode' : mode,
                    'text': text,
                    }
                cls.db.append(litmus)

    @staticmethod
    def count():
        return len(Litmus.db)

    @staticmethod
    def get_by_id(id):
        litmus = Litmus.db[id]
        return litmus

    @staticmethod
    def classify_by_group(sort, order):
        db = {}
        for star in Litmus.supernova:
            db.update({star['litmus']['name']: {'count':0, 'litmus':[] }})
        for litmus in Litmus.db:
            group = litmus['group']
            db[group]['count'] = db[group]['count'] + 1
            db[group]['litmus'].append(litmus)
        keys = db.keys()
        for group in keys:
            if order == "ascend":
                db[group]['litmus'] = sorted(db[group]['litmus'], key=lambda g: g[sort])
            elif order == "descend":
                db[group]['litmus'] = sorted(db[group]['litmus'], reverse=True, key=lambda g: g[sort])
        
        return db

    @staticmethod
    def get_text(mode) :
        if mode == 'Light' :
            text_color = '#000000'
            text_font = 'bold'
        else :
            text_color = '#FFFFFF'
            text_font = 'normal'
        return {'color': text_color, 'font': text_font}

    
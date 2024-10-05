# BioSeq_control

## О BioSeq_control
***Bioseq_control*** пресдтавляет для небольшую программу для парсинга и манипуляции сиквенсов.
В программу входят два модуля и главный скрипт.


#### Структура файлов:


- **main_script.py**
- **modules_HW4**
    - **dna_rna_tools**
    - **fastq_filter_m**

### Установка

```bash
git clone git@github.com:goomtreker/BioSEQ_control.git
```
Подключение модуля:
Запустите любой интерпритатор pyton
```pyton
import sys
sys.path.append("path/to/BioSEQ_control")
```
### Использование
В главном скрипте ```main_script``` лежит две функции:
-```run_dna_rna_tools```
-```filter_fastq```

```run_dna_rna_tool``` принимает одну или несколько сиквенсов (РНК или ДНК)
 и последним примнимается агрумент, который будет идентичным названию
 функции.
возможные функции:
- ```reverse```
- ```reverse_complement```
- ```transcribe```
- ```complement```
- ```find_possible_ORF```
- ```GC_status```
- ```check_acid_type```
Пример использования:
```python
run_dna_rna_tools('GGCccttggATC', 'gcccggttt', 'reverse_complement')
```

```filter_fastq``` на вход принимает словарь вида: {Имя сиквенса: (сиквенс, символы качества рида)}. Пример:
```
EXAMPLE_FASTQ = {
    '@SRX079801': ('ACAGCAACATAAACATGA....AA', 'FGGGFGGGFGGGFG...GGD'),
    '@SRX079802': ('ATTA...TG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    ...
    '@SRX079812': ('AGT...AATGACCCG', '<98;...BBC')
    }      
```
А также пороги по длине, гц составу, и контролю качетва рида.
На выходе будут сиквенсы, которые прошли все пороги.
Пример использования:
```python
filter_fastq(
        seqs: EXAMPLE_FASTQ,
        gc_bounds=44,
        length_bounds=(10000, 100000), quality_threshold=34):
```

Примечание:
Вопросы замечания и предложения можно задать/предложить во вкладке issues.
P.S Нас часто спрашивают знаем ли мы что-нибудь о пакете Bio.SeqIO?
Наш ответ:
Мы никогда о нем не слышали.:)
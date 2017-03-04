#ifndef Q_ROOT_CANVAS_H
#define Q_ROOT_CANVAS_H

#include <QWidget>
#include "TCanvas.h"

class TObject;
class QObject;
class QEvent;
class QPaintEvent;
class QResizeEvent;
class QMouseEvent;

class QRootCanvas : public QWidget
{
    Q_OBJECT

public:
    QRootCanvas(QWidget *parent = 0);
    virtual ~QRootCanvas();
    void Refresh() {fCanvas->Modified(); fCanvas->Update();};
    TCanvas *GetCanvas() {return fCanvas;};

    // wrapper class
    void SetFillColor(Color_t c);
    void SetFrameFillColor(Color_t c);
    void SetGrid();
    void cd(Int_t subpadnumber=0);
    virtual void Browse(TBrowser *b);
    void Clear(Option_t *option="");
    void Close(Option_t *option="");
    virtual void Draw(Option_t *option="");
    virtual TObject *DrawClone(Option_t *option="");
    virtual TObject *DrawClonePad();
    virtual void EditorBar();
    void EnterLeave(TPad *prevSelPad, TObject *prevSelObj);
    void FeedbackMode(Bool_t set);
    void Flush();
    void UseCurrentStyle();
    void ForceUpdate();
    const char *GetDISPLAY();
    TContextMenu *GetContextMenu();
    Int_t GetDoubleBuffer();
    Int_t GetEvent();
    Int_t GetEventX();
    Int_t GetEventY();
    Color_t GetHighLightColor();
    TVirtualPad *GetPadSave();
    TObject *GetSelected();
    Option_t *GetSelectedOpt();
    TVirtualPad *GetSelectedPad();
    Bool_t GetShowEventStatus();
    Bool_t GetAutoExec();
    Size_t GetXsizeUser();
    Size_t GetYsizeUser();
    Size_t GetXsizeReal();
    Size_t GetYsizeReal();
    Int_t GetCanvasID();
    Int_t GetWindowTopX();
    Int_t GetWindowTopY();
    UInt_t GetWindowWidth();
    UInt_t GetWindowHeight();
    UInt_t GetWw();
    UInt_t GetWh();
    virtual void GetCanvasPar(Int_t &wtopx, Int_t &wtopy, UInt_t &ww, UInt_t &wh);
    virtual void HandleInput(EEventType button, Int_t x, Int_t y);
    Bool_t HasMenuBar();
    void Iconify();
    Bool_t IsBatch();
    Bool_t IsRetained();
    virtual void ls(Option_t *option="");
    void MoveOpaque(Int_t set=1);
    Bool_t OpaqueMoving();
    Bool_t OpaqueResizing();
    virtual void Paint(Option_t *option="");
    virtual TPad *Pick(Int_t px, Int_t py, TObjLink *&pickobj);
    virtual TPad *Pick(Int_t px, Int_t py, TObject *prevSelObj);
    virtual void Resize(Option_t *option="");
    void ResizeOpaque(Int_t set=1);
    void SaveSource(const char *filename="", Option_t *option="");
    virtual void SetCursor(ECursor cursor);
    virtual void SetDoubleBuffer(Int_t mode=1);
    void SetWindowPosition(Int_t x, Int_t y);
    void SetWindowSize(UInt_t ww, UInt_t wh);
    void SetCanvasSize(UInt_t ww, UInt_t wh);
    void SetHighLightColor(Color_t col);
    void SetSelected(TObject *obj);
    void SetSelectedPad(TPad *pad);
    void Show();
    virtual void Size(Float_t xsizeuser=0, Float_t ysizeuser=0);
    void SetBatch(Bool_t batch=kTRUE);
    void SetRetained(Bool_t retained=kTRUE);
    void SetTitle(const char *title="");
    virtual void ToggleEventStatus();
    virtual void ToggleAutoExec();
    virtual void Update();

signals:
    void TObjectSelected(TObject *, TCanvas *);

protected:
    TCanvas *fCanvas;

    virtual bool eventFilter(QObject *o, QEvent *e);
    virtual void resizeEvent(QResizeEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);
    virtual void mousePressEvent(QMouseEvent *e);
    virtual void mouseReleaseEvent(QMouseEvent *e);
    virtual void mouseDoubleClickEvent(QMouseEvent *e);
    virtual void paintEvent(QPaintEvent *e);
    virtual void leaveEvent(QEvent *e);

private:
    bool fNeedResize;
};

#endif

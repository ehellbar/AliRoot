// $Id$
// Category: interfaces
//
// Author: D. Adamova
//
//========================================================
//
//------------TG4MainFrame.cxx--------------------------------//
//---------Main Window for the AG4 Geometry Browser---//
//
//========================================================= 

#include "TG4MainFrame.h"
#include "TG4Editor.h"
#include "TG4ListTreeFrame.h"
#include "TG4VolumesFrames.h"
#include "TG4MaterialsFrames.h"
#include "TG4GuiVolume.h"
#include "TG4Globals.h"

#include <TGListTree.h>
#include <TGTab.h>
#include <TGMenu.h>
#include <TApplication.h>
#include <TGMsgBox.h>
#include <TGTextBuffer.h>

#include <G4LogicalVolume.hh>

ClassImp(TG4MainFrame)

TG4MainFrame::TG4MainFrame(const TGWindow* p, UInt_t w, UInt_t h)
    : TGMainFrame(p, w, h)
{
//---> Creates the main frame 

  fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
  fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

  fPopupMenu= new TGPopupMenu(gClient->GetRoot());
  fPopupMenu->AddEntry("&Close Window", 1);
  fPopupMenu->AddEntry("&Exit Root", 2);
  
  fPopupMenuTest = new TGPopupMenu(this);
  fPopupMenuTest->AddEntry("&Message", 3);
    
  fPopupMenuHelp = new TGPopupMenu(gClient->GetRoot());
  fPopupMenuHelp->AddEntry("&About", 4);

   fPopupMenu->Associate(this);
   fPopupMenuTest->Associate(this);
   fPopupMenuHelp->Associate(this);
 
   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&CloseWindow/ExitRoot",fPopupMenu, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Draw Control", fPopupMenuTest, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Report", fPopupMenuHelp, fMenuBarHelpLayout);

   AddFrame(fMenuBar, fMenuBarLayout);
 
//------>Adding a tab
   fTab = new TGTab(this, 400, 400);
   TGLayoutHints* lTabLayout = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1);
   AddFrame(fTab, lTabLayout);

//------->Frame for ListTree of logical volumes
   TGCompositeFrame* tf = fTab->AddTab("Volumes");
   flistTreeFrame = new TG4ListTreeFrame( tf, this);
   fVolumesListTree=flistTreeFrame->GetVolumesListTree();

//----->Frame for volumes properties
   tf = fTab->AddTab("Volumes Properties");
   fvolumesFrames = new TG4VolumesFrames( tf, this);


//----->Frame for materials properties
   tf = fTab->AddTab("Materials Properties");
   fmaterialsFrames = new TG4MaterialsFrames( tf, this);

//----->Window name and final mapping
   SetWindowName("ALICE Geant4 Browser");
   MapSubwindows();
   Resize(GetDefaultSize());
 
   MapWindow();
   
}

TG4MainFrame::TG4MainFrame(const TG4MainFrame& mf)
   : TGMainFrame( (const TGMainFrame&) mf)
{
// Dummy copy constructor 
  TG4Globals::Exception(
    "Attempt to use TG4MainFrame copy constructor.");
}

TG4MainFrame& TG4MainFrame::operator=(const TG4MainFrame& mf)
{
  // check assignement to self
  if (this == &mf) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4MainFrame singleton.");
    
  return *this;  
}

TG4MainFrame::~TG4MainFrame()
{
//----> Delete created widgets.
   G4cout << "\n Now in the TG4MainFrame destructor\n" << G4endl;  
   delete fMenuBarLayout;
   delete fMenuBarItemLayout;
   delete fMenuBarHelpLayout;

   delete fPopupMenu;
   delete fPopupMenuTest;
   delete fPopupMenuHelp;
 
   delete fMenuBar;
   delete fTab;

   delete fVolumesListTree;
   delete flistTreeFrame;
   delete fvolumesFrames;   
   delete fmaterialsFrames;
}

TG4VolumesFrames* TG4MainFrame::GetVolumesFrames() const
{
//---> For use in TG4GeometryGUI
   return fvolumesFrames;
}

TG4MaterialsFrames* TG4MainFrame::GetMaterialsFrames() const
{
//---> For use in TG4GeometryGUI
   return fmaterialsFrames;
}

TGListTreeItem*  TG4MainFrame::
AddItem(TObject* obj, TGListTreeItem* parent, const char* name, 
                       const TGPicture* open, const TGPicture* closed)
{
//----->Add item to the list tree
    return fVolumesListTree->AddItem(parent, name, obj, open, closed);
}

void TG4MainFrame::CloseWindow()
{
   // Got close message for this MainFrame. Calls parent CloseWindow()
   // (which destroys the window) ((//and terminate the application)).
   // The close message is generated by the window manager when its close
   // window menu item is selected.

   TGMainFrame::CloseWindow();
   gApplication->Terminate(0);
}

Bool_t TG4MainFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
//---> Processes events generated by the widgets in the frame

//=============================================================
//-----> text buffers for Message Boxes
   TGTextBuffer* lMsgBTtleBf = new TGTextBuffer(100);
   TGTextBuffer* lMsgBAnnBf1 = new TGTextBuffer(100);
   TGTextBuffer* lMsgBAnnBf2 = new TGTextBuffer(100);
   lMsgBTtleBf->AddText(0, "MsgBox");
   lMsgBAnnBf1->AddText(0, "Volumes drawing almost completed !");
   lMsgBAnnBf2->AddText(0, "YOU'RE CLOSING THE MAIN WINDOW!");


//=================================================================
//----->Message Box

   int  retval;
   EMsgBoxIcon icontype = kMBIconExclamation;
   Int_t buttons;
//===================================================================
//----->Editor window text
const char *editortxt =
"This is a to-be AG4 Geometry Browser. \n"
"The volumes drawing coming soon.\n" ;


//=================================================================
//----->Process messages to widgets
    switch (GET_MSG(msg)) {
    
    case kC_TEXTENTRY:
        switch (GET_SUBMSG(msg)) {
        case kTE_TEXTCHANGED:
            switch (parm1) {
            case 300:
                G4cout <<" Acting in TextEntry !!! " << G4endl;
		
                break;
		
             default:
	         break;
            };
	    
	 default:
	    break;
         }
        break;

//---->case Handle Popup menus    
    case kC_COMMAND:
        switch (GET_SUBMSG(msg)) {
            case kCM_MENU:
               switch (parm1) {

                   case 1:
		   
		     buttons = kMBOk;
                    // for (Int_t i=1; i<3; i++)
                    //      buttons |= i;
		     new TGMsgBox(fClient->GetRoot(), this,
                                  lMsgBTtleBf->GetString(), lMsgBAnnBf2->GetString(),
                                  icontype, buttons, &retval);
		           // if not here, produces
		           // Error in <RootX11ErrorHandler>: BadWindow 
		           // (invalid Window parameter) (XID: 100663461)		      
                     TGMainFrame::CloseWindow();
                             break;  

                  case 2:
		     G4cout << "\n\n!!!!EXITING. BYE!!!\n\n" << G4endl; 
                     CloseWindow();   //-->and exit root
                     break;
		     
		   case 3:
		      //-->popup Message 
		      buttons = kMBDismiss;
		      new TGMsgBox(fClient->GetRoot(), this,
                                  lMsgBTtleBf->GetString(), lMsgBAnnBf1->GetString(),
                                  icontype, buttons, &retval);
		      break;
		      
		  case 4:
		     //-->editor window
		     {
		      TG4Editor* ed = new TG4Editor(this, 400, 150);
                      ed->LoadBuffer(editortxt);
		      ed->Popup();
		      };
		      break;
		  
		  default:    
		     break;
                };
		
            case kCM_COMBOBOX: 
	        switch(parm1) {
		
		    case 100:
		      fvolumesFrames->DisplayVolumeCharacteristics();
		      break;
		      
		    case 200:
		      fmaterialsFrames->DisplayMaterialCharacteristics();
		      break;
		     
                    default:
                      break;
               };
	       
            default:
               break;
         }
	 break;	   

//----->case Handle volumes ListTree
    case kC_LISTTREE:
	switch (GET_SUBMSG(msg)) {

//----->Cases to Handle mouse click
   //-->case 1 
	case kCT_ITEMCLICK: 
    //---> Button 1: Select volume
            if (parm1 == kButton1){
	    TGListTreeItem* item = fVolumesListTree->GetSelected();
             if (item) {
	     
	     TG4GuiVolume* volume=((TG4GuiVolume*) item->GetUserData());   
	     G4LogicalVolume* lvolume = volume->GetLogicalVolume();

	      if  ( lvolume ) {
	       G4cout << "The selected logical volume name is   " 
	            << lvolume->GetName() << G4endl;
		  };
	      };
	    }; 
	      
    //---> Button 3: Draw Volume
            if (parm1 == kButton3){
             TGListTreeItem* item = fVolumesListTree->GetSelected();
	     if (item){
                flistTreeFrame->DrawSelectedVolume(item);};
	    }; 
		  
	    break;

   //-->case 2	    
	case kCT_ITEMDBLCLICK:
	    if (parm1 == kButton1) {
		if (fVolumesListTree->GetSelected() != 0) {
		    gClient->NeedRedraw(fVolumesListTree);
		};
	    };
	    break;
   //-->default for GET_SUBMSG	    
	default:
	    break;
	}
	break;
//---->default for GET_MSG	
    default:
	break;
    }

    delete lMsgBTtleBf;
    delete lMsgBAnnBf1;
    delete lMsgBAnnBf2; 
   
    return kTRUE;
}

